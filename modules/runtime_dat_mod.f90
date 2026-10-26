! =====================================================================
! runtime_dat_mod.f90 —— runtime.dat v3 完整消费端 (工程主通路)
!   读入: SPECIES/IC/MET(7)/BC/PHOT(l,m,n)/EMIS/DRY/WET/DIL/SCEN/SERIES
!   求值: 气象插值(含 BLH 导数) · SZA→J (MCM 参数化) · 排放日型 ·
!         情景因子(排放组缩放/JFAC/DT) · M/O2/H2O 环境
!   设计: 纯数据模块 (无 CVODE 依赖), 逐行 token 解析, 供多驱动共享
! =====================================================================
module runtime_dat_mod
    use oz_mech, only: ns, nj
    implicit none
    integer, parameter :: NMET = 7        ! TEMP PRES RH RADIAT BLH RAIN JFAC
    integer, parameter :: NSCEN_MAX = 16, NSER_MAX = 8, NTOK = 260
    real(8), parameter :: PI = 3.14159265358979d0

    character(len=16) :: sp_name(ns)
    real(8) :: y0(ns)
    real(8) :: days, dtout, rtol_cfg, atol_cfg
    real(8) :: met_c(NMET) = (/298.15d0, 1013.25d0, 50.0d0, 0.0d0, 1000.0d0, 0.0d0, 1.0d0/)
    integer :: met_s(NMET) = 0
    real(8) :: tau_mix = 86400.0d0, bc_c(ns)
    integer :: bc_s(ns)
    real(8) :: j_l(nj), j_m(nj), j_n(nj)
    real(8) :: emis0(ns), emis_amp(ns)
    integer :: emis_box(ns)               ! 1=BOX 日型 0=FLAT
    real(8) :: dry_v(ns)
    integer :: dry_on(ns)
    real(8) :: wet_a(ns), wet_b(ns)
    real(8) :: dil_k = 0.0d0
    real(8) :: phot_lat = 45.0d0, phot_doy = 172.0d0
    integer :: phot_mode = 1     ! R16: 1=MCM(l·cos^m·exp(-n secχ)) 2=DIRECT(l 恒定, 不调制)
    integer :: nscen = 0, nser = 0
    character(len=32) :: scen_name(NSCEN_MAX)
    real(8) :: scen_fnox(NSCEN_MAX) = 1.0d0, scen_fvoc(NSCEN_MAX) = 1.0d0, &
               scen_fco(NSCEN_MAX) = 1.0d0, scen_fch4(NSCEN_MAX) = 1.0d0, &
               scen_jfac(NSCEN_MAX) = 1.0d0, scen_dt(NSCEN_MAX) = 0.0d0
    integer :: ser_n(NSER_MAX) = 0
    real(8) :: ser_t(NSER_MAX, 128), ser_v(NSER_MAX, 128)
    ! 情景生效态 (OpenMP 线程私有: 每线程可独立设情景)
    real(8) :: emis_scale(ns) = 1.0d0, jfac_scen = 1.0d0, dt_scen = 0.0d0
    !$omp threadprivate(emis_scale, jfac_scen, dt_scen)
    ! 情景分组 (单源化: 原 batch_mod 双实现收敛于此; SO2 归 VOC 组—含硫排放
    ! 与氧化性 VOC 同源谱, schema 不增组的取舍)
    integer :: grpnox(ns) = 0, grpvoc(ns) = 0, grpco(ns) = 0, grpch4(ns) = 0

contains
    ! ---------- 行 token 切分 ----------
    subroutine split(line, tok, n)
        character(len=*), intent(in) :: line
        character(len=32), intent(out) :: tok(NTOK)
        integer, intent(out) :: n
        character(len=len(line)) :: s
        integer :: p
        s = adjustl(line)
        n = 0
        do while (len_trim(s) > 0 .and. n < NTOK)
            p = index(s, " ")
            if (p == 0) p = len_trim(s) + 1
            n = n + 1
            tok(n) = s(1:p - 1)
            s = adjustl(s(p:))
        end do
    end subroutine split

    subroutine read_runtime(path)
        character(len=*), intent(in) :: path
        character(len=4096) :: line   ! V2.3: 74 物种 BC 行 ~1117 字符 > 旧 1024
        character(len=32) :: tok(NTOK), sec
        integer :: u, ios, nt, i, k, idx, isc, iser, nscen_expect

        k = 0
        ios = 0; sec = ""; isc = 0; iser = 0; nscen_expect = 0
        open (newunit=u, file=path, status="old", action="read")
        do while (ios == 0)
            read (u, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(adjustl(line)) == 0) cycle
            if (index(adjustl(line), "#") == 1) cycle
            call split(line, tok, nt)
            if (nt == 0) cycle
            select case (sec)
            case ("")
                ! ---- 段头 ----
                select case (tok(1))
                case ("RUN")
                    ! R11-P3: 校验 ns 与 SPECIES 行一致 (旧实现读过就丢)
                    if (int_val(tok(2)) /= ns) then
                        print '(A,I0,A,I0,A)', "[read_runtime FATAL] RUN.ns=", &
                            int_val(tok(2)), " != 生成码 ns=", ns, &
                            " (机制/生成码/runtime 三者不同步)"
                        error stop
                    end if
                    nscen_expect = int_val(tok(5))   ! RUN ns nj nser nscen
                case ("SPECIES")
                    do i = 1, min(ns, nt - 1)
                        sp_name(i) = trim(tok(i + 1))
                    end do
                case ("MODEL")
                    days = num(tok(2)); dtout = num(tok(3))
                    rtol_cfg = num(tok(4)); atol_cfg = num(tok(5))
                case ("DIL")
                    sec = "DIL"
                case ("PHOT")
                    ! PHOT 头行: MODE lat doy (后续 nj 行 J idx l m n)
                    ! 行只有 3 个 token: MODE=1, lat=2, doy=3 → 判据必须是 nt>=3
                    ! (旧写 nt>=4 使 DOY 永远读不到, 静默停留在默认值 172)
                    if (nt >= 3) phot_lat = num(tok(2))
                    if (nt >= 3) phot_doy = num(tok(3))
                    sec = "PHOT"
                case ("IC")
                    sec = "IC"
                case ("MET")
                    sec = "MET"; k = 0
                case ("BC")
                    sec = "BC"
                case ("EMIS")
                    sec = "EMIS"; k = 0
                case ("DRY")
                    sec = "DRY"; k = 0
                case ("WET")
                    sec = "WET"; k = 0
                case ("SCEN")
                    sec = "SCEN"
                case ("SERIES")
                    sec = "SERIES"
                end select
            case ("IC")
                do i = 1, min(ns, nt)
                    y0(i) = num(tok(i))
                end do
                sec = ""
            case ("MET")
                k = k + 1
                if (tok(1) == "C") then
                    met_c(k) = num(tok(2)); met_s(k) = 0
                else if (tok(1) == "S") then
                    met_s(k) = int_val(tok(2))
                end if
                if (k == NMET) sec = ""
            case ("BC")
                if (nt /= 1 + 2*ns) then
                    print *, "[read_runtime FATAL] BC 行 token 数 ", nt, " != ", 1 + 2*ns, " (行截断?)"
                    error stop
                end if
                tau_mix = num(tok(1))
                do i = 1, ns
                    bc_c(i) = num(tok(1 + i))
                    bc_s(i) = int_val(tok(1 + ns + i))
                end do
                sec = ""
            case ("PHOT")
                if (tok(1) == "MCM" .or. tok(1) == "DIRECT") then
                    phot_mode = merge(1, 2, tok(1) == "MCM")   ! R16: 记录模式
                    if (nt >= 3) phot_lat = num(tok(2))
                    if (nt >= 3) phot_doy = num(tok(3))
                else if (tok(1) == "J") then
                    idx = int_val(tok(2))
                    j_l(idx) = num(tok(3)); j_m(idx) = num(tok(4)); j_n(idx) = num(tok(5))
                    if (idx == nj) sec = ""
                end if
            case ("EMIS")
                k = k + 1
                emis0(k) = num(tok(2))
                emis_box(k) = merge(1, 0, tok(3) == "BOX")
                emis_amp(k) = num(tok(4))
                if (k == ns) sec = ""
            case ("DRY")
                k = k + 1
                dry_v(k) = num(tok(2)); dry_on(k) = int_val(tok(3))
                if (k == ns) sec = ""
            case ("WET")
                k = k + 1
                wet_a(k) = num(tok(2)); wet_b(k) = num(tok(3))
                if (k == ns) sec = ""
            case ("DIL")
                dil_k = num(tok(1))
                sec = ""
            case ("SCEN")
                isc = isc + 1
                if (isc <= NSCEN_MAX) then
                    scen_name(isc) = trim(tok(1))
                    scen_fnox(isc) = num(tok(2)); scen_fvoc(isc) = num(tok(3))
                    scen_fco(isc) = num(tok(4)); scen_fch4(isc) = num(tok(5))
                    scen_jfac(isc) = num(tok(6)); scen_dt(isc) = num(tok(7))
                end if
                if (isc == nscen_expect) sec = ""
            case ("SERIES")
                ! 行1: name n; 行2: t...; 行3: v...
                iser = iser + 1
                if (iser > NSER_MAX) then
                    print '(A,I0,A)', "[read_runtime FATAL] 序列数 ", iser, &
                        " > NSER_MAX=8 (cfgbuild 应已拦截; runtime.dat 被手改?)"
                    error stop
                end if
                nser = iser
                ser_n(iser) = int_val(tok(2))
                if (ser_n(iser) > 128) then
                    print '(A,I0,A)', "[read_runtime FATAL] 序列点数 ", ser_n(iser), &
                        " > 128 (ser_t 上限; cfgbuild 应已拦截)"
                    error stop
                end if
                if (trim(tok(1)) == "BLH") met_s(5) = iser
                if (trim(tok(1)) == "RAIN") met_s(6) = iser
                read (u, '(A)', iostat=ios) line
                if (ios /= 0) exit
                call split(line, tok, nt)
                do i = 1, min(ser_n(iser), nt)
                    ser_t(iser, i) = num(tok(i))
                end do
                read (u, '(A)', iostat=ios) line
                if (ios /= 0) exit
                call split(line, tok, nt)
                do i = 1, min(ser_n(iser), nt)
                    ser_v(iser, i) = num(tok(i))
                end do
            end select
        end do
        close (u)
        nscen = max(isc, 0)
        if (nscen == 0) then
            nscen = 1; scen_name(1) = "DEFAULT"
        end if
        ! R11-P3: 缺 SPECIES 行 → 硬失败 (旧实现兜底为 SPi, 致情景组映射静默全失效)
        do i = 1, ns
            if (len_trim(sp_name(i)) == 0) then
                print '(A)', "[read_runtime FATAL] runtime.dat 缺 SPECIES 行 (物种名不可用," &
                    // " 情景组缩放会静默失效) — 请重跑 cfgbuild"
                error stop
            end if
        end do
    contains
        real(8) function num(s)
            character(len=*), intent(in) :: s
            read (s, *) num
        end function
        integer function int_val(s)
            character(len=*), intent(in) :: s
            read (s, *) int_val
        end function
    end subroutine read_runtime

    ! ---------- 序列线性插值 (外延取端值) ----------
    pure function ser_interp(k, t) result(v)
        integer, intent(in) :: k
        real(8), intent(in) :: t
        real(8) :: v
        integer :: i
        v = ser_v(k, 1)
        do i = 1, ser_n(k) - 1
            if (t <= ser_t(k, i + 1)) then
                v = ser_v(k, i) + (ser_v(k, i + 1) - ser_v(k, i))* &
                    (t - ser_t(k, i))/max(ser_t(k, i + 1) - ser_t(k, i), 1.0d0)
                return
            end if
        end do
        v = ser_v(k, ser_n(k))
    end function ser_interp

    ! ---------- 气象求值 (BLH 单位 m→cm; dBLH 为 cm/s) ----------
    subroutine eval_met(t, TEMP, PRES, RH, BLH, RAIN, JFAC, dBLH)
        real(8), intent(in) :: t
        real(8), intent(out) :: TEMP, PRES, RH, BLH, RAIN, JFAC, dBLH
        TEMP = val_or_ser(1, t) + dt_scen
        PRES = val_or_ser(2, t)
        RH = val_or_ser(3, t)
        BLH = val_or_ser(5, t)*1.0d2
        RAIN = max(val_or_ser(6, t), 0.0d0)
        JFAC = val_or_ser(7, t)
        dBLH = (val_or_ser(5, t + 600.0d0) - val_or_ser(5, t - 600.0d0))*1.0d2/1200.0d0
    contains
        real(8) function val_or_ser(k, tt)
            integer, intent(in) :: k
            real(8), intent(in) :: tt
            if (met_s(k) > 0 .and. met_s(k) <= nser) then
                val_or_ser = ser_interp(met_s(k), tt)
            else
                val_or_ser = met_c(k)
            end if
        end function
    end subroutine eval_met

    ! ---------- 环境 M / O2 / H2O ----------
    subroutine eval_env(TEMP, PRES, RH, M, O2, H2O)
        real(8), intent(in) :: TEMP, PRES, RH
        real(8), intent(out) :: M, O2, H2O
        real(8) :: es, e
        M = 6.022d23*(PRES*1.0d2)/(8.314d0*TEMP)*1.0d-6
        O2 = 0.2095d0*M
        es = 6.112d0*exp(17.67d0*(TEMP - 273.15d0)/(TEMP - 29.65d0))
        e = (RH/100.0d0)*es
        H2O = 6.022d23*(e*1.0d2)/(8.314d0*TEMP)*1.0d-6
    end subroutine eval_env

    ! ---------- SZA → J: J = jf·l·cos^m(χ)·exp(−n·sec(χ)) ----------
    subroutine eval_J(t, JARR)
        real(8), intent(in) :: t
        real(8), intent(out) :: JARR(nj)
        real(8) :: cosz, decl, lat, ha, hour, jf
        integer :: i
        hour = mod(t, 86400.0d0)/3600.0d0
        decl = 23.45d0*sin(2.0d0*PI*(phot_doy - 81.0d0)/365.0d0)*PI/180.0d0
        lat = phot_lat*PI/180.0d0
        ha = 2.0d0*PI*(hour - 12.0d0)/24.0d0
        cosz = sin(lat)*sin(decl) + cos(lat)*cos(decl)*cos(ha)
        if (met_s(7) > 0 .and. met_s(7) <= nser) then
            jf = jfac_scen*ser_interp(met_s(7), t)
        else
            jf = jfac_scen*met_c(7)
        end if
        do i = 1, nj
            if (cosz > 1.0d-3 .and. j_l(i) > 0.0d0) then
                if (phot_mode == 1) then
                    ! MCM: 天顶角调制 l·cos^m·exp(-n/cosχ)
                    JARR(i) = jf*j_l(i)*cosz**j_m(i)*exp(-j_n(i)/cosz)
                else
                    ! DIRECT (R16): 不施加天顶角调制 —— 白天恒定取 l, 夜间为 0
                    JARR(i) = jf*j_l(i)
                end if
            else
                JARR(i) = 0.0d0
            end if
        end do
    end subroutine eval_J

    ! ---------- 排放日型 BOX: 0.1+0.9·amp·sin²(π(h−6)/12) ----------
    pure function box_shape(t, amp)
        real(8), intent(in) :: t, amp
        real(8) :: box_shape, hh
        hh = mod(t, 86400.0d0)
        if (hh >= 6.0d0*3600.0d0 .and. hh <= 18.0d0*3600.0d0) then
            box_shape = 0.1d0 + 0.9d0*amp* &
                sin(PI*(hh - 6.0d0*3600.0d0)/(12.0d0*3600.0d0))**2
        else
            box_shape = 0.1d0
        end if
    end function box_shape

    ! ---------- 情景生效 ----------
    subroutine build_groups()
        ! 排放→情景组分映射 (单源; 新类 V2.2 全入 VOC 组)
        integer :: i
        grpnox = 0; grpvoc = 0; grpco = 0; grpch4 = 0
        do i = 1, ns
            select case (trim(sp_name(i)))
            case ("NO", "NO2")
                grpnox(i) = 1
            case ("HCHO", "CH3CHO", "ISOP", "DMS", "MVK", "MACR", "MSA", "SO2", &
                  "ALK", "OLE", "ARO1", "ARO2", "KET", "TERP")
                grpvoc(i) = 1
            case ("CO")
                grpco(i) = 1
            case ("CH4")
                grpch4(i) = 1
            end select
        end do
    end subroutine build_groups

    subroutine scen_build(isc)
        ! R10-P0: 单源组映射 —— 只读 build_groups 的 grp* 数组, 不再维护第二份物种清单
        ! (旧实现漏缩 SO2/ALK/OLE/ARO1/ARO2/KET/TERP, 致 fac_v3 情景沙盘对主排放源失效;
        !  batch 侧 esc_of 早已用 grp*, 两条路线由此同源)
        integer, intent(in) :: isc
        integer :: i
        call build_groups()
        emis_scale = 1.0d0
        jfac_scen = scen_jfac(isc)
        dt_scen = scen_dt(isc)
        do i = 1, ns
            if (grpnox(i) == 1) emis_scale(i) = scen_fnox(isc)
            if (grpvoc(i) == 1) emis_scale(i) = scen_fvoc(isc)
            if (grpco(i)  == 1) emis_scale(i) = scen_fco(isc)
            if (grpch4(i) == 1) emis_scale(i) = scen_fch4(isc)
        end do
    end subroutine scen_build
end module runtime_dat_mod
