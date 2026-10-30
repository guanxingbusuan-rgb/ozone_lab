! =====================================================================
! fac_batch_omp.f90 —— FAC-V3 批量响应面 · OpenMP 并行 · 指示因子全家族
!   扫描: (f_NOx × f_VOC × jfac × ΔT) 排放/气象网格 (F2: 后两维 CLI 可选)
!   诊断: 末 24 步（要求 post-spinup ≥24 步）逐时均值 + 12–18h 午后口径
!         (F5, *_pm 列, Sillman 原文口径；批量为聚合输出，不落全时间序列)
!   ⚠ 统计窗 = **末 24 个输出步**（dtout=3600 时 = 模拟最后 1 天）：
!     O3max/O3mean/HCHOmax/NO2mean 均为该窗内统计，**不是 8h 日最大**。
!   官方指示因子 5 族: FNR, H2O2/HNO3, L_NOx/Q (Kleinman 2005),
!                      P(H2O2)/P(HNO3), OPE   + 状态量 NOy
!   备用扩展 (仍占列, 不进判据): O3/NOy, HCHO/NOy, O3/NOz
!   IPR (F1): 每案例末 24h 窗口均值 → <outf>_ipr.csv (fnox,fvoc,jfac,dt+ns×6)
!   状态侧车 (实验 J): <outf>_state.csv —— 同窗口内 NO2/HCHO 的**一致均值对**
!     (全窗 + 12–18h 双口径)，供状态空间 EKMA (x=NO2, y=HCHO, 等值线=O3) 使用。
!     主 CSV 的 HCHOmax/NO2mean 一个取极值一个取均值, 不能当两轴用。
!   预算表: batch_budget_init 启动期解析索引 (自检硬失败); 氮守恒在线不变量
!   物理层: chem_phys_mod 唯一实现 (A1, 与 fac_v3 共享; F3 序列 BC 在此接线)
!   线程模型: C API 直连, 每线程独立 mem/N_Vector/user_data
! CLI: fac_batch_omp [ng≤64] [out.csv] [runtime.dat] [jfac列表] [dt列表] [fmin] [fmax]
!      [模式 EMIS|IC] [第9参] [LIN]
!      缺省网格 fmin=0.3、步长比 1.5 (等比)；传 fmin/fmax 指定范围，第 10 参 LIN 走线性网格
!      正式入口 = EKMA-FNR / EKMA-emission 的 run.py（传绝对路径）；手工等价例:
!        fac_batch_omp 11 raw.csv runtime_fnr.dat 1.0 0.0 1 11 IC NOX:HCHO LIN
!        fac_batch_omp 26 raw.csv runtime_emission.dat 1.0 0.0 1 7.25 EMIS 1 LIN
! =====================================================================
module batch_mod
    use iso_c_binding, only: c_double, c_int
    use oz_mech, only: ns, nr, nq, nj, mre, mpe, update_rates
    use runtime_dat_mod
    use chem_phys_mod, only: phys_tend, eval_met_dt, eval_J_jf, NPHYS
    ! 情景分组用 runtime_dat_mod 的 grpnox/grpvoc/grpco/grpch4 (单源化)
    implicit none
    integer :: nreact(nr), reac_id(mre, nr)
    integer :: nprod(nr), prod_id(mpe, nr)
    real(8) :: prod_coef(mpe, nr)
    integer :: nnoy = 0
    integer :: noy_id(32) = 0, noy_w(32) = 0   ! R10-P2: 容量 20→32
    integer, parameter :: NT_Q = 6, NT_LN = 2
    integer :: q_rx(NT_Q) = 0, ln_rx(NT_LN) = 0
    real(8) :: q_cf(NT_Q) = 0.0d0
    integer :: npo3 = 0, nphno3 = 0, nph2o2 = 0, npnoz = 0
    integer :: po3_rx(mpe*nr), phno3_rx(mpe*nr), ph2o2_rx(mpe*nr), pnoz_rx(mpe*nr)
    real(8) :: po3_cf(mpe*nr), phno3_cf(mpe*nr), ph2o2_cf(mpe*nr), pnoz_cf(mpe*nr)
    integer :: spO3 = 0, spHCHO = 0, spNO2 = 0, spH2O2 = 0, spHNO3 = 0, spNO = 0
    real(8) :: wn(ns) = 0.0d0
    integer :: ic1_mask(ns) = 0, ic2_mask(ns) = 0 ! IC 扰动物种掩膜 (a9 = ic1spec:ic2spec)
    logical :: ic_mode = .false.                  ! IC 扰动模式开关
    logical :: ndbg_once = .true.
    !$omp threadprivate(ndbg_once)
    integer :: nNviol = 0                     ! 氮守恒不变量违例案例数（共享，critical 内累加）
    real(8) :: peff_last(nr)               ! diag_indicators 最近一次 peff (IRR 累积源)
    !$omp threadprivate(peff_last)
    character(len=1), parameter :: iprl(6) = ["P", "L", "E", "D", "B", "X"]

    type, bind(C) :: runp_t
        real(c_double) :: fnox = 1.0d0, fvoc = 1.0d0, fco = 1.0d0, fch4 = 1.0d0
        real(c_double) :: jfac = 1.0d0, dts = 0.0d0
        ! IC 扰动模式 (实验 J --ic): 两个 IC 组的缩放因子 (O3 等因变量不缩放)
        real(c_double) :: ic1 = 1.0d0, ic2 = 1.0d0
        integer(c_int) :: icmode = 0
    end type
contains
    integer function spidx_(nm)
        character(len=*), intent(in) :: nm
        integer :: i
        spidx_ = 0
        do i = 1, ns
            if (trim(sp_name(i)) == trim(nm)) then
                spidx_ = i; return
            end if
        end do
    end function spidx_

    subroutine batch_init(g)
        character(len=*), intent(in) :: g
        integer :: u, k, i
        open (newunit=u, file=g//"/oz.reac", status="old")
        do k = 1, nr
            read (u, *) nreact(k), (reac_id(i, k), i=1, nreact(k))
        end do
        close (u)
        open (newunit=u, file=g//"/oz.prod", status="old")
        do k = 1, nr
            read (u, *) nprod(k), (prod_id(i, k), prod_coef(i, k), i=1, nprod(k))
        end do
        close (u)
        call build_groups()
        call noy_add("NO", 1); call noy_add("NO2", 1); call noy_add("HNO3", 1)
        call noy_add("HONO", 1); call noy_add("NO3", 1); call noy_add("N2O5", 2)
        call noy_add("PAN", 1); call noy_add("MPAN", 1)
        call noy_add("ISOPNIT", 1); call noy_add("ISOPDNT", 1)
        call noy_add("ClONO2", 1); call noy_add("BrONO2", 1)
        call noy_add("ISOPNO3", 1)
        call noy_add("PNA", 1)      ! V2.3: HO2NO2 温度计
        call noy_add("RONO2", 1)    ! V2.3: ALK 烷基硝酸酯 (15% 分支)
        call noy_add("OLENIT", 1)   ! V2.3: OLE 硝酸酯 (3% 分支)
        call noy_add("TERPNIT", 1)  ! V2.3: 单萜硝酸酯 (NO3 加成)
        ! V2.4: 新增加成/硝酸酯家族 (芳香-甲酚/烷烃/醚/酯/倍半萜/乙烯/内烯烃/PPAN)
        call noy_add("CSLNIT", 1)
        call noy_add("LALKNIT", 1)
        call noy_add("MTBENIT", 1)
        call noy_add("ESTERNIT", 1)
        call noy_add("SESQNIT", 1)
        call noy_add("ETHNIT", 1)
        call noy_add("IOLENIT", 1)
        call noy_add("PPAN", 1)
    end subroutine

    subroutine noy_add(nm, w)
        character(len=*), intent(in) :: nm
        integer, intent(in) :: w
        integer :: i
        do i = 1, ns
            if (trim(sp_name(i)) == trim(nm)) then
                if (nnoy >= size(noy_id)) then
                    print '(A,A,A,I0)', "[FATAL] NOy 物种数超上限: ", trim(nm), &
                        " (noy_id 容量 ", size(noy_id), ")"
                    error stop
                end if
                nnoy = nnoy + 1
                noy_id(nnoy) = i; noy_w(nnoy) = w
                return
            end if
        end do
        print '(A,A)', "  [警告] NOy 物种不在机制: ", trim(nm)
    end subroutine

    ! ---------- 预算表 (启动期一次解析) ----------
    subroutine batch_budget_init()
        integer :: idxs(nr), n, k, i, j, isp, m
        call collect_rxn(["O1D   "], ["OH    ", "OH    "], idxs, n)
        if (n > 0) then; q_rx(1) = idxs(1); q_cf(1) = 2.0d0; end if
        call collect_rxn(["HONO  "], ["NO    ", "OH    "], idxs, n)
        if (n > 0) then; q_rx(2) = idxs(1); q_cf(2) = 1.0d0; end if
        call collect_rxn(["H2O2  "], ["OH    ", "OH    "], idxs, n)
        if (n > 0) then; q_rx(3) = idxs(1); q_cf(3) = 2.0d0; end if
        call collect_rxn(["HCHO  "], ["HO2   ", "HO2   ", "CO    "], idxs, n)
        if (n > 0) then; q_rx(4) = idxs(1); q_cf(4) = 2.0d0; end if
        call collect_rxn(["CH3OOH"], ["CH3O  ", "OH    "], idxs, n)
        if (n > 0) then; q_rx(5) = idxs(1); q_cf(5) = 1.0d0; end if
        call collect_rxn(["ISOPOOH"], ["ISOPBO", "OH    "], idxs, n)
        if (n > 0) then; q_rx(6) = idxs(1); q_cf(6) = 1.0d0; end if
        call collect_rxn(["NO2   ", "OH    "], ["HNO3  "], idxs, n)
        if (n > 0) ln_rx(1) = idxs(1)
        call collect_rxn(["NO    ", "OH    "], ["HONO  "], idxs, n)
        if (n > 0) ln_rx(2) = idxs(1)
        npo3 = build_prod("O3", po3_rx, po3_cf)
        nphno3 = build_prod("HNO3", phno3_rx, phno3_cf)
        nph2o2 = build_prod("H2O2", ph2o2_rx, ph2o2_cf)
        m = 0
        do j = 1, nnoy
            select case (trim(sp_name(noy_id(j))))
            case ("NO", "NO2", "N2O5")
                cycle
            end select
            isp = noy_id(j)
            do k = 1, nr
                do i = 1, nprod(k)
                    if (prod_id(i, k) == isp) then
                        m = m + 1
                        pnoz_rx(m) = k; pnoz_cf(m) = prod_coef(i, k)
                    end if
                end do
            end do
        end do
        npnoz = m
        spO3 = spidx_("O3"); spHCHO = spidx_("HCHO"); spNO2 = spidx_("NO2")
        spH2O2 = spidx_("H2O2"); spHNO3 = spidx_("HNO3"); spNO = spidx_("NO")
        ! R10-P2: 索引 0 = 未找到 → yp(0) 越界读, 必须硬失败
        if (min(spO3, spHCHO, spNO2, spH2O2, spHNO3, spNO) <= 0) then
            print '(A,A,A)', "[FATAL] 机制缺关键物种: O3=", trim(sp_name(spO3)), &
                " (批次诊断需要 O3/HCHO/NO2/H2O2/HNO3/NO 全部存在)"
            error stop
        end if
        wn = 0.0d0
        do j = 1, nnoy
            wn(noy_id(j)) = real(noy_w(j), 8)
        end do
    end subroutine batch_budget_init

    integer function build_prod(nm, rx, cf) result(nn)
        character(len=*), intent(in) :: nm
        integer, intent(out) :: rx(mpe*nr)
        real(8), intent(out) :: cf(mpe*nr)
        integer :: k, i, isp
        isp = spidx_(nm)
        nn = 0
        if (isp <= 0) return
        do k = 1, nr
            do i = 1, nprod(k)
                if (prod_id(i, k) == isp) then
                    nn = nn + 1
                    rx(nn) = k; cf(nn) = prod_coef(i, k)
                end if
            end do
        end do
    end function build_prod

    ! ---------- 匹配核心 ----------
    logical function rxn_match(k, lhs, rhs)
        integer, intent(in) :: k
        character(len=*), intent(in) :: lhs(:), rhs(:)
        integer :: i, c1, c2
        rxn_match = .false.
        if (nreact(k) /= size(lhs) .or. nprod(k) /= size(rhs)) return
        do i = 1, size(lhs)
            call tally(lhs(i), c1, c2)
            if (c1 /= cnt_reac(k, lhs(i)) .or. c2 /= cnt_prod(k, lhs(i))) return
        end do
        do i = 1, size(rhs)
            call tally(rhs(i), c1, c2)
            if (c1 /= cnt_reac(k, rhs(i)) .or. c2 /= cnt_prod(k, rhs(i))) return
        end do
        rxn_match = .true.
    contains
        subroutine tally(nm, c1, c2)
            character(len=*), intent(in) :: nm
            integer, intent(out) :: c1, c2
            c1 = count(lhs == nm)
            c2 = count(rhs == nm)
        end subroutine
        integer function cnt_reac(k, nm) result(c)
            integer, intent(in) :: k
            character(len=*), intent(in) :: nm
            integer :: i
            c = 0
            do i = 1, nreact(k)
                if (trim(sp_name(reac_id(i, k))) == trim(nm)) c = c + 1
            end do
        end function
        integer function cnt_prod(k, nm) result(c)
            integer, intent(in) :: k
            character(len=*), intent(in) :: nm
            integer :: i
            c = 0
            do i = 1, nprod(k)
                if (trim(sp_name(prod_id(i, k))) == trim(nm)) c = c + 1
            end do
        end function
    end function rxn_match

    subroutine collect_rxn(lhs, rhs, idxs, n)
        character(len=*), intent(in) :: lhs(:), rhs(:)
        integer, intent(out) :: idxs(nr), n
        integer :: k
        n = 0
        do k = 1, nr
            if (rxn_match(k, lhs, rhs)) then
                n = n + 1
                idxs(n) = k
            end if
        end do
    end subroutine collect_rxn

    ! ---------- 速率核 + 物理倾向 (A1: 物理唯一实现在 chem_phys_mod) ----------
    subroutine batch_rates(t, y, rp, peff)
        real(8), intent(in) :: t, y(ns)
        type(runp_t), intent(in) :: rp
        real(8), intent(out) :: peff(nr)
        real(8) :: q(nq), p(nr), JARR(nj)
        real(8) :: TEMP, PRES, RH, BLH, RAIN, JF, dB, MM, O2, H2O
        integer :: k, i
        call eval_met_dt(t, rp%dts, TEMP, PRES, RH, BLH, RAIN, JF, dB)
        call eval_env(TEMP, PRES, RH, MM, O2, H2O)
        JF = rp%jfac*JF          ! R11-P0: 绝对因子 = CLI 网格 jfac × config JFAC
        call eval_J_jf(t, JF, JARR)   ! eval_J_jf 内部会归一化掉 config 侧, 不再二次相乘
        call update_rates(TEMP, MM, O2, H2O, JARR, y, q, p)
        peff = p
        do k = 1, nr
            do i = 1, nreact(k)
                peff(k) = peff(k)*y(reac_id(i, k))
            end do
        end do
    end subroutine batch_rates

    subroutine batch_resid(t, y, ydot, rp)
        real(8), intent(in) :: t, y(ns)
        real(8), intent(out) :: ydot(ns)
        type(runp_t), intent(in) :: rp
        real(8) :: peff(nr), yphys(ns), iphys(ns, NPHYS), metv(7)
        integer :: k, i
        call batch_rates(t, y, rp, peff)
        call phys_tend(t, y, esc_of(rp), yphys, iphys, metv)
        ydot = yphys
        do k = 1, nr
            do i = 1, nreact(k)
                ydot(reac_id(i, k)) = ydot(reac_id(i, k)) - peff(k)
            end do
            do i = 1, nprod(k)
                ydot(prod_id(i, k)) = ydot(prod_id(i, k)) + prod_coef(i, k)*peff(k)
            end do
        end do
    end subroutine batch_resid

    ! ---------- IC 扰动模式 (实验 J --ic / EKMA-FNR): 精确物种对 ----------
    !   a9 = "ic1spec:ic2spec"，每个 spec 为:
    !     物种名    → 只缩放该物种的初值 (如 "NO2"、"HCHO")
    !     NOX       → NO + NO2 两组 (NOx 总量语义)
    !     VOC       → 活性 VOC 组 (HCHO/CH3CHO/ISOP/MVK/MACR/MGLY)
    !   O3 永不缩放 (因变量); 未列物种保持基线。
    subroutine ic_spec_mask(spec, msk)
        character(len=*), intent(in) :: spec
        integer, intent(out) :: msk(ns)
        integer :: i, n
        msk = 0
        select case (trim(spec))
        case ("NOX", "nox")
            do i = 1, ns
                if (trim(sp_name(i)) == "NO" .or. trim(sp_name(i)) == "NO2") msk(i) = 1
            end do
        case ("VOC", "voc")
            do i = 1, ns
                select case (trim(sp_name(i)))
                case ("HCHO", "CH3CHO", "ISOP", "MVK", "MACR", "MGLY"); msk(i) = 1
                end select
            end do
        case default
            do i = 1, ns
                if (trim(sp_name(i)) == trim(spec)) msk(i) = 1
            end do
            n = sum(msk)
            if (n == 0) then
                print '(3A)', "[FATAL] IC 扰动 spec 未匹配任何机制物种: '", trim(spec), "'"
                error stop
            end if
        end select
    end subroutine ic_spec_mask

    subroutine apply_ic_scale(y, rp)
        real(8), intent(inout) :: y(ns)
        type(runp_t), intent(in) :: rp
        integer :: k
        if (rp%icmode /= 1) return
        do k = 1, ns
            if (ic1_mask(k) == 1) then
                y(k) = y(k)*rp%ic1
            else if (ic2_mask(k) == 1) then
                y(k) = y(k)*rp%ic2
            end if
        end do
    end subroutine apply_ic_scale

    ! 情景缩放 → esc 数组 (A2: 网格语义)
    function esc_of(rp) result(esc)
        type(runp_t), intent(in) :: rp
        real(8) :: esc(ns)
        integer :: k
        esc = 1.0d0
        if (rp%icmode == 1) then
            esc = 0.0d0          ! IC 扰动模式 = 纯初值化学 (经典 EKMA 约定)
            return
        end if
        do k = 1, ns
            if (grpnox(k) == 1) esc(k) = rp%fnox
            if (grpvoc(k) == 1) esc(k) = rp%fvoc
            if (grpco(k) == 1) esc(k) = rp%fco
            if (grpch4(k) == 1) esc(k) = rp%fch4
        end do
    end function esc_of

    ! ---------- 指示因子 + IPR 累积 (pm=.true. 时仅 12–18h, F5) ----------
    subroutine diag_indicators(t, y, rp, d, iprc, pm, acc_ipr)
        real(8), intent(in) :: t, y(ns)
        type(runp_t), intent(in) :: rp
        real(8), intent(out) :: d(15)
        real(8), intent(inout) :: iprc(ns, 6)
        logical, intent(in) :: pm
        logical, intent(in), optional :: acc_ipr
        logical :: do_ipr
        real(8) :: peff(nr), noy, noz, qrad, lnox, ph2o2, phno3, po3, pnoz
        real(8) :: ych(ns), yphys(ns), iphys(ns, NPHYS), metv(7), nsum, l1
        real(8) :: hh
        integer :: i, k
        d = 0.0d0
        do_ipr = .true.
        if (present(acc_ipr)) do_ipr = acc_ipr
        if (pm) then
            hh = mod(t, 86400.0d0)/3600.0d0
            if (hh < 12.0d0 .or. hh > 18.0d0) return   ! 窗外: d=0 (intent(out) 安全)
        end if
        call batch_rates(t, y, rp, peff)
        peff_last = peff
        ! 氮守恒在线不变量 (每案例首快照; ndbg_once 在案例开头重置, threadprivate ⇒ 无竞态)
        if (ndbg_once) then
            ndbg_once = .false.
            ych = 0.0d0
            do k = 1, nr
                do i = 1, nreact(k)
                    ych(reac_id(i, k)) = ych(reac_id(i, k)) - peff(k)
                end do
                do i = 1, nprod(k)
                    ych(prod_id(i, k)) = ych(prod_id(i, k)) + prod_coef(i, k)*peff(k)
                end do
            end do
            nsum = dot_product(wn, ych)
            l1 = sum(wn*abs(ych))
            if (abs(nsum) > 1.0d-4*max(l1, 1.0d0)) then
                !$omp critical (nwarn)
                nNviol = nNviol + 1
                print '(A,ES10.3)', "  [氮守恒 FAIL] |Σ w_N·dy/dt| = ", abs(nsum)
                !$omp end critical (nwarn)
            end if
        end if
        ! IPR 累积: 机制 P/L + 物理四账 (R10-P1: 仅首次调用累积, 防 12-18h 窗双计)
        if (do_ipr) then
        do k = 1, nr
            do i = 1, nreact(k)
                iprc(reac_id(i, k), 2) = iprc(reac_id(i, k), 2) - peff(k)
            end do
            do i = 1, nprod(k)
                iprc(prod_id(i, k), 1) = iprc(prod_id(i, k), 1) + prod_coef(i, k)*peff(k)
            end do
        end do
        call phys_tend(t, y, esc_of(rp), yphys, iphys, metv)
        do k = 1, ns
            iprc(k, 3) = iprc(k, 3) + iphys(k, 1)
            iprc(k, 4) = iprc(k, 4) + iphys(k, 2)
            iprc(k, 5) = iprc(k, 5) + iphys(k, 3)
            iprc(k, 6) = iprc(k, 6) + iphys(k, 4)
        end do
        end if   ! do_ipr
        ! 浓度比 (第一类)
        noy = 0.0d0
        do i = 1, nnoy
            noy = noy + noy_w(i)*max(y(noy_id(i)), 0.0d0)
        end do
        noz = noy - max(y(spNO), 0.0d0) - max(y(spNO2), 0.0d0)
        ! 过程量 (第二类, 预算表热路径)
        qrad = 0.0d0
        do i = 1, NT_Q
            if (q_rx(i) > 0) qrad = qrad + q_cf(i)*peff(q_rx(i))
        end do
        lnox = 0.0d0
        do i = 1, NT_LN
            if (ln_rx(i) > 0) lnox = lnox + peff(ln_rx(i))
        end do
        ph2o2 = 0.0d0
        do i = 1, nph2o2
            ph2o2 = ph2o2 + ph2o2_cf(i)*peff(ph2o2_rx(i))
        end do
        phno3 = 0.0d0
        do i = 1, nphno3
            phno3 = phno3 + phno3_cf(i)*peff(phno3_rx(i))
        end do
        po3 = 0.0d0
        do i = 1, npo3
            po3 = po3 + po3_cf(i)*peff(po3_rx(i))
        end do
        pnoz = 0.0d0
        do i = 1, npnoz
            pnoz = pnoz + pnoz_cf(i)*peff(pnoz_rx(i))
        end do
        d(1) = rat(y(spH2O2), y(spHNO3))
        d(2) = rat(y(spO3), noy)
        d(3) = rat(y(spHCHO), noy)
        d(4) = rat(y(spO3), noz)
        d(5) = noy/2.46d10
        d(6) = rat(y(spHCHO), y(spNO2))
        d(13) = lnox; d(14) = ph2o2; d(15) = phno3
        d(10) = qrad; d(11) = po3; d(12) = pnoz
    contains
        real(8) function rat(a, b)
            real(8), intent(in) :: a, b
            rat = 0.0d0
            if (b > 1.0d-30) rat = a/b
        end function
    end subroutine diag_indicators
end module

! ---- C 回调 ----
integer(c_int) function cv_rhs(t, yv, ydv, ud) bind(C)
    use batch_mod
    use cvode_c_mod, only: c_nvgetarr
    use iso_c_binding
    real(c_double), value :: t
    type(c_ptr), value :: yv, ydv, ud
    real(8), pointer :: y(:), yd(:)
    type(runp_t), pointer :: rp
    call c_f_pointer(c_nvgetarr(yv), y, [ns])
    call c_f_pointer(c_nvgetarr(ydv), yd, [ns])
    call c_f_pointer(ud, rp)
    call batch_resid(t, y, yd, rp)
    cv_rhs = 0_c_int
end function

program fac_batch_omp
    use batch_mod
    use cvode_c_mod
    use iso_c_binding
    use omp_lib
    implicit none
    procedure(cvrhs_proc) :: cv_rhs
    integer :: ng, i, j, k, nsteps, nspin, u, uipr, u3, ust, ier, cnt, cntpm, ijm, njf, ndt, it, ip1
    integer :: iw, isp, ir, ilen             ! 侧车写出循环变量 (勿复用 i —— 曾致 IRR 坐标全错)
    character(len=256) :: a1, outf, rtfile, a4, a5, a8, a9, a10, a6, a7
    character(len=16) :: numbuf                ! 单个数字的临时格式 (线程私有)
    character(len=12000) :: bufi               ! IPR 行缓冲 (线程私有, 临界区外格式化)
    character(len=4096) :: bufr                ! IRR 行缓冲
    real(8) :: tstart, tend, fs(64), fjs(16), dts(16), TRET, fmin, fmax, gr
    type(c_ptr), target :: mem, ny, nat
    type(runp_t), target :: rp
    real(8), pointer :: yp(:), ap(:)
    integer(c_long) :: nst, nst_acc
    integer :: nfail = 0, nwarn = 0
    real(8) :: o3max, hchomax, no2m, o3m, dsum(15), dg(15)
    ! 状态空间 EKMA 用的一致统计量 (同一窗口内各物种/过程量各自均值)
    real(8) :: hchom, no2pm, hchopm, hno3m, h2o2m, hno3pm, h2o2pm
    real(8) :: dpms(15)
    real(8) :: ipracc(ns, 6)
    real(8) :: irracc(nr)
    real(8) :: O3v, HCHOv, NO2v, HNO3v, H2O2v, hh
    integer :: iO3v, iHCHOv, iNO2v, iHNO3v, iH2O2v

    call get_command_argument(1, a1); call get_command_argument(2, outf)
    call get_command_argument(3, rtfile); call get_command_argument(4, a4)
    call get_command_argument(5, a5)
    call get_command_argument(6, a6); call get_command_argument(7, a7)
    call get_command_argument(8, a8); call get_command_argument(9, a9)
    call get_command_argument(10, a10)
    if (len_trim(a1) == 0) a1 = "9"
    if (len_trim(outf) == 0) outf = "../output/batch_v3.csv"
    if (len_trim(rtfile) == 0) rtfile = "../runtime_v3.dat"
    if (len_trim(a4) == 0) a4 = "1.0"
    if (len_trim(a5) == 0) a5 = "0.0"
    read (a1, *) ng
    ! 网格几何: 缺省保持旧网格 fmin=0.3, 步长比 1.5 (哨兵/文档兼容);
    ! 传第 6/7 参数可指定 f 范围, 加密同范围分辨率 (实验 J: 13 点覆盖 0.3–7.7)
    fmin = 0.3d0
    fmax = 0.3d0*1.5d0**(ng - 1)
    if (len_trim(a6) > 0) read (a6, *) fmin
    if (len_trim(a7) > 0) read (a7, *) fmax
    if (ng > size(fs)) then
        print '(A,I0,A,I0,A)', "[FATAL] ng=", ng, " 超过网格数组上限 ", size(fs)
        error stop
    end if
    if (trim(a10) == "LIN" .or. trim(a10) == "lin") then
        ! 线性网格: f = fmin + (fmax-fmin)*(i-1)/(ng-1)
        ! EKMA 正交扰动实验用 (0%→200% 以 5% 等距 ⇒ f=1.00→3.00, 41 点)
        do i = 1, ng
            fs(i) = fmin
            if (ng > 1) fs(i) = fmin + (fmax - fmin)*dble(i - 1)/dble(ng - 1)
            if (fs(i) <= 0.0d0) then
                print '(A,F9.4,A)', "[FATAL] 线性网格出现非正因子 ", fs(i), " (fmin 必须 > 0)"
                error stop
            end if
        end do
    else
        gr = 1.0d0
        if (ng > 1) gr = (fmax/fmin)**(1.0d0/(ng - 1))
        do i = 1, ng
            fs(i) = fmin*gr**(i - 1)
        end do
    end if
    njf = split_list(a4, fjs)
    ndt = split_list(a5, dts)
    ! 第 8 参 = 模式 (EMIS 缺省 | IC); 第 9 参 = IC 第二组 (1=VOC 2=ISOP 3=CO 4=CH4)
    ic_mode = (trim(a8) == "IC" .or. trim(a8) == "ic")   ! a9 = "ic1spec:ic2spec"

    call read_runtime(trim(rtfile))
    call batch_init("gen_v3")
    call batch_budget_init()
    if (ic_mode) then
        ip1 = index(a9, ":")
        if (ip1 == 0) then
            print '(A)', "[FATAL] IC 模式第 9 参须为 'ic1spec:ic2spec' (如 NO2:HCHO / NOX:VOC)"
            error stop
        end if
        call ic_spec_mask(a9(1:ip1-1), ic1_mask)
        call ic_spec_mask(a9(ip1+1:), ic2_mask)
        tau_mix = 1.0d30        ! IC 模式关闭 BC 张弛 (纯初值化学)
        print '(5A,I0,A,I0,A)', "  [IC 模式] 缩放 ", trim(a9(1:ip1-1)), " × ", &
            trim(a9(ip1+1:)), "; 排放关闭; BC 张弛关闭; 网格 ", ng, "×", ng
    end if
    iO3v = spO3; iHCHOv = spHCHO; iNO2v = spNO2
    iHNO3v = spHNO3; iH2O2v = spH2O2
    nsteps = nint(days*86400.0d0/dtout)
    nspin = nint(spinup_cfg*86400.0d0/dtout)
    if (nsteps - nspin < 24) then
        print '(A,F5.2,A,F5.2,A)', "FAC-BATCH FATAL: post-spinup only ", &
            (nsteps - nspin)*dtout/86400.0d0, " d (SPINUP=", spinup_cfg, &
            " d), need >= 24 output steps for the 24-step aggregate"
        stop 1
    end if
    call diag_selfcheck()

    call cpu_time(tstart)
    open (newunit=u, file=trim(outf), status="replace")
    write (u, '(A)') "fnox,fvoc,jfac,dt,O3max,O3mean,HCHOmax,NO2mean,FNR_mean,"// &
        "H2O2HNO3,O3NOy,HCHONOy,O3NOz,NOy_ppb,LnOxQ,PH2O2PHNO3,OPE,"// &
        "FNR_pm,H2O2HNO3_pm,O3NOy_pm,HCHONOy_pm,O3NOz_pm,NOy_pm,nst"
    open (newunit=uipr, file=trim(outf)//"_ipr.csv", status="replace")
    open (newunit=u3, file=trim(outf)//"_irr.csv", status="replace")   ! F1-batch: 窗口均值 IRR
    write (u3, '(A)', advance="no") "fnox,fvoc,jfac,dt"
    do k = 1, nr
        write (u3, '(",r", I0)', advance="no") k
    end do
    write (u3, *)
    write (uipr, '(A)', advance="no") "fnox,fvoc,jfac,dt"
    do k = 1, ns
        do i = 1, 6
            write (uipr, '(",", A, "_", A)', advance="no") trim(sp_name(k)), iprl(i)
        end do
    end do
    write (uipr, *)
    ! 状态空间 EKMA 侧车 (实验 J): 与主 CSV 同案例顺序, 每案例一行一致统计量。
    ! 通用状态表: 物种均值 ppb + 过程量均值 molec/cm3/s, 全窗与 12–18h 双口径,
    ! 供 --axes 任意拆分子分母构造响应面 (S2–S5)。
    open (newunit=ust, file=trim(outf)//"_state.csv", status="replace")
    write (ust, '(A)') "fnox,fvoc,jfac,dt,"// &
        "NO2mean_ppb,HCHOmean_ppb,HNO3mean_ppb,H2O2mean_ppb,"// &
        "NO2mean_pm_ppb,HCHOmean_pm_ppb,HNO3mean_pm_ppb,H2O2mean_pm_ppb,"// &
        "Q,L_NOx,P_H2O2,P_HNO3,P_O3,P_NOz,"// &
        "Q_pm,L_NOx_pm,P_H2O2_pm,P_HNO3_pm,P_O3_pm,P_NOz_pm,"// &
        "O3max_ppb,O3mean_ppb,nst"

    !$omp parallel do private(ijm, i, j, k, it, mem, ny, nat, rp, yp, ap, TRET, nst, &
    !$omp& nst_acc, bufi, bufr, numbuf, ilen, iw, isp, ir, &
    !$omp& ier, cnt, cntpm, o3max, o3m, hchomax, hchom, no2m, no2pm, hchopm, &
    !$omp& hno3m, h2o2m, hno3pm, h2o2pm, dsum, dpms, dg, &
    !$omp& ipracc, irracc, O3v, HCHOv, NO2v, HNO3v, H2O2v, hh) &
    !$omp& reduction(+:nfail, nwarn) schedule(dynamic) collapse(2)
    do it = 1, njf*ndt
    do ijm = 1, ng*ng
        i = (ijm - 1)/ng + 1
        j = mod(ijm - 1, ng) + 1
        rp%fnox = fs(i); rp%fvoc = fs(j)
        rp%fco = 1.0d0; rp%fch4 = 1.0d0
        rp%ic1 = fs(i); rp%ic2 = fs(j)
        rp%icmode = merge(1_c_int, 0_c_int, ic_mode)
        rp%jfac = fjs((it - 1)/ndt + 1); rp%dts = dts(mod(it - 1, ndt) + 1)
        ny = c_nvnew(int(ns, c_long))
        nat = c_nvnew(int(ns, c_long))
        call c_f_pointer(c_nvgetarr(ny), yp, [ns])
        call c_f_pointer(c_nvgetarr(nat), ap, [ns])
        yp = y0
        call apply_ic_scale(yp, rp)
        ap = atol_cfg
        call tune_atol_(ap)
        mem = c_cvcreate(CV_BDF, CV_NEWTON)
        ier = c_cvinit(mem, c_funloc(cv_rhs), 0.0d0, ny)
        ier = c_cvsvtol(mem, rtol_cfg, nat)
        ier = c_cvdense(mem, int(ns, c_long))
        ier = c_cvsetmaxsteps(mem, 200000_c_long)
        ier = c_cvsetmaxstep(mem, 100.0d0)
        ier = c_cvsetud(mem, c_loc(rp))
        nst_acc = 0_c_long
        o3max = 0.0d0; o3m = 0.0d0; hchomax = 0.0d0; no2m = 0.0d0
        hchom = 0.0d0; no2pm = 0.0d0; hchopm = 0.0d0
        hno3m = 0.0d0; h2o2m = 0.0d0; hno3pm = 0.0d0; h2o2pm = 0.0d0
        dsum = 0.0d0; dpms = 0.0d0
        ipracc = 0.0d0
        irracc = 0.0d0
        cnt = 0; cntpm = 0
        ndbg_once = .true.    ! 每案例都检一次首快照（threadprivate ⇒ 线程内私有，无竞态；
                              ! 若只在每线程检一次，fatal 化会变成"看哪个线程先碰哪个案例"）
        do k = 1, nsteps
            ier = c_cvode(mem, dble(k)*dtout, ny, TRET, CV_NORMAL)
            if (ier /= 0) then
                ! 韧性重启 (与 fac_v3 同策略): 同容差 + MAX_STEP 10 s,
                ! 从当前状态重新初始化; 极端角点 (低 NOx 高 VOC) 收敛失败可救回
                ier = c_cvgetnumsteps(mem, nst)
                nst_acc = nst_acc + nst
                call c_cvfree(c_loc(mem))
                mem = c_cvcreate(CV_BDF, CV_NEWTON)
                ier = c_cvinit(mem, c_funloc(cv_rhs), max(TRET, 0.0d0), ny)
                ier = c_cvsvtol(mem, rtol_cfg, nat)
                ier = c_cvdense(mem, int(ns, c_long))
                ier = c_cvsetmaxsteps(mem, 200000_c_long)
                ier = c_cvsetmaxstep(mem, 10.0d0)
                ier = c_cvsetud(mem, c_loc(rp))
                ier = c_cvode(mem, max(dble(k)*dtout, TRET + 1.0d0), ny, TRET, CV_NORMAL)
                if (ier /= 0) then
                    ! 已积累统计的案例: 降级为警告 (极端角点收敛失败, 行仍写出)
                    if (cnt > 0) then
                        ! 已有统计的失败: 行仍写出但窗口不完整 → 计入降级, 供调用方记录
                        nwarn = nwarn + 1
                        !$omp critical (warn_print)
                        print '(2(A,F7.3),A,I0,A)', "  [警告] 案例 fnox=", rp%fnox, &
                            " fvoc=", rp%fvoc, " 在 t=", k, " 步收敛失败 (保留已积累统计)"
                        !$omp end critical (warn_print)
                    else
                        nfail = nfail + 1
                    end if
                    exit
                end if
            end if
            if (k > nsteps - 24) then
                O3v = yp(iO3v)/2.46d10; HCHOv = yp(iHCHOv)/2.46d10; NO2v = yp(iNO2v)/2.46d10
                HNO3v = yp(iHNO3v)/2.46d10; H2O2v = yp(iH2O2v)/2.46d10
                o3max = max(o3max, O3v); o3m = o3m + O3v
                hchomax = max(hchomax, HCHOv); no2m = no2m + NO2v
                hchom = hchom + HCHOv
                hno3m = hno3m + HNO3v; h2o2m = h2o2m + H2O2v
                call diag_indicators(dble(k)*dtout, yp, rp, dg, ipracc, .false.)
                dsum = dsum + dg
                irracc(1:nr) = irracc(1:nr) + peff_last(1:nr)
                cnt = cnt + 1
                hh = mod(dble(k)*dtout, 86400.0d0)/3600.0d0
                if (hh >= 12.0d0 .and. hh <= 18.0d0) then
                    ! 同一状态: 午后口径的诊断量与全窗完全相同 → 复用 dg, 不再二次算速率
                    dpms = dpms + dg
                    no2pm = no2pm + NO2v; hchopm = hchopm + HCHOv
                    hno3pm = hno3pm + HNO3v; h2o2pm = h2o2pm + H2O2v
                    cntpm = cntpm + 1
                end if
            end if
        end do
        ier = c_cvgetnumsteps(mem, nst)
        nst = nst_acc + nst
        call c_cvfree(c_loc(mem))
        call c_nvdestroy(ny); call c_nvdestroy(nat)
        if (cnt > 0) then
            ! 侧车行在线程内格式化（临界区外），临界区只做整行写 —— 降低写盘串行化
            write (bufi, '(3(F6.3, ","), F6.1)') fs(i), fs(j), rp%jfac, rp%dts
            ilen = len_trim(bufi)
            do isp = 1, ns
                do iw = 1, 6
                    write (numbuf, '(ES15.5E3)') ipracc(isp, iw)/cnt
                    bufi(ilen + 1:ilen + 1 + len_trim(numbuf)) = ","//trim(numbuf)
                    ilen = ilen + 1 + len_trim(numbuf)
                end do
            end do
            write (bufr, '(3(F6.3, ","), F6.1)') fs(i), fs(j), rp%jfac, rp%dts
            ilen = len_trim(bufr)
            do ir = 1, nr
                write (numbuf, '(ES15.5E3)') irracc(ir)/cnt
                bufr(ilen + 1:ilen + 1 + len_trim(numbuf)) = ","//trim(numbuf)
                ilen = ilen + 1 + len_trim(numbuf)
            end do
            !$omp critical (csv_write)
            write (u, '(4(F6.3, ","), 4(F7.2, ","), 9(F7.3, ","), 6(F8.3, ","), I0)') &
                fs(i), fs(j), rp%jfac, rp%dts, &
                o3max, o3m/cnt, hchomax, no2m/cnt, &
                dsum(6)/cnt, dsum(1)/cnt, dsum(2)/cnt, dsum(3)/cnt, &
                dsum(4)/cnt, dsum(5)/cnt, rat0(dsum(13), dsum(10)), &
                rat0(dsum(14), dsum(15)), rat0(dsum(11), dsum(12)), &
                dpms(6)/max(cntpm, 1), dpms(1)/max(cntpm, 1), dpms(2)/max(cntpm, 1), &
                dpms(3)/max(cntpm, 1), dpms(4)/max(cntpm, 1), dpms(5)/max(cntpm, 1), nst
            write (ust, '(4(F6.3, ","), 8(F9.4, ","), 6(ES15.5E3, ","), ' // &
                '6(ES15.5E3, ","), 2(F9.4, ","), I0)') &
                fs(i), fs(j), rp%jfac, rp%dts, &
                no2m/cnt, hchom/cnt, hno3m/cnt, h2o2m/cnt, &
                no2pm/max(cntpm, 1), hchopm/max(cntpm, 1), &
                hno3pm/max(cntpm, 1), h2o2pm/max(cntpm, 1), &
                dsum(10)/cnt, dsum(13)/cnt, dsum(14)/cnt, dsum(15)/cnt, &
                dsum(11)/cnt, dsum(12)/cnt, &
                dpms(10)/max(cntpm, 1), dpms(13)/max(cntpm, 1), &
                dpms(14)/max(cntpm, 1), dpms(15)/max(cntpm, 1), &
                dpms(11)/max(cntpm, 1), dpms(12)/max(cntpm, 1), &
                o3max, o3m/cnt, nst
            write (uipr, '(A)') trim(bufi)
            write (u3, '(A)') trim(bufr)
            !$omp end critical (csv_write)
        end if
    end do
    end do
    !$omp end parallel do
    close (u); close (uipr); close (u3); close (ust)
    call cpu_time(tend)
    print '(A,I0,A,I0,A,I0,A,I0,A,F10.2,A,I0,A)', "FAC-BATCH-OMP ok (", ng*ng*njf*ndt, &
        " 案例, 失败 ", nfail, ", 降级 ", nwarn, ", 氮守恒违例 ", nNviol, ", ", &
        tend - tstart, " CPU s, 线程 ", omp_get_max_threads(), ")"
    ! 求解失败或化学不自洽都必须让调用方看见：截断网格 / NOy 不守恒的响应面都不能出
    if (nfail > 0 .or. nNviol > 0) then
        print '(A,I0,A,I0,A)', "[FATAL] 求解失败 ", nfail, " 例 / 氮守恒违例 ", nNviol, &
            " 例 —— 结果不完整或化学不自洽，调用方应拒绝出表"
        stop 1
    end if
contains
    integer function split_list(s, arr) result(n)
        character(len=*), intent(in) :: s
        real(8), intent(out) :: arr(16)
        character(len=32) :: tok(32)
        integer :: nt, i
        call splitcsv(s, tok, nt)
        n = 0
        do i = 1, min(nt, 16)
            read (tok(i), *) arr(i)
            n = n + 1
        end do
    end function
    subroutine splitcsv(s, tok, n)
        character(len=*), intent(in) :: s
        character(len=32), intent(out) :: tok(32)
        integer, intent(out) :: n
        integer :: p
        character(len=len(s)) :: w
        w = s
        n = 0
        do
            p = index(w, ",")
            if (p == 0) then
                n = n + 1; tok(n) = adjustl(w); exit
            end if
            n = n + 1; tok(n) = adjustl(w(1:p - 1))
            w = w(p + 1:)
        end do
    end subroutine
    real(8) function rat0(a, b)
        real(8), intent(in) :: a, b
        rat0 = 0.0d0
        if (b > 1.0d-30) rat0 = a/b
    end function
    subroutine diag_selfcheck()
        integer :: nmiss
        nmiss = 0
        if (npo3 == 0) nmiss = nmiss + 1
        if (nphno3 == 0) nmiss = nmiss + 1
        if (nph2o2 == 0) nmiss = nmiss + 1
        if (npnoz == 0) nmiss = nmiss + 1
        if (any(q_rx == 0)) nmiss = nmiss + 1
        if (any(ln_rx == 0)) nmiss = nmiss + 1
        if (nmiss > 0) then
            print '(A,I0,A)', "  [警告] ", nmiss, " 组预算项未解析 — 检查机制名!"
            stop 1
        else
            print '(A,I0,A,I0,A,I0,A,I0,A)', "  [自检] 预算表解析 ✓ Q×", count(q_rx > 0), &
                " + L_NOx×", count(ln_rx > 0), " + 产率 ", npo3 + nphno3 + nph2o2 + npnoz, &
                " 项 (NOy ", nnoy, " 物种)"
        end if
    end subroutine
    subroutine tune_atol_(a)
        ! 单源转调 chem_phys_mod (V2.3 起与 fac_v3 共表, 含 TERP/RONO2 等新物种)
        use chem_phys_mod, only: tune_atol_by_name
        real(8), intent(inout) :: a(ns)
        call tune_atol_by_name(a)
    end subroutine
end program
