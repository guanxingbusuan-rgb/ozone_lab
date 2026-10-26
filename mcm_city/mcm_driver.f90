! =====================================================================
! mcm_driver.f90 —— 真实 MCM v3.3.1 CH₄ 子集 · 城市光化学 5 天模拟
!   FAC-SPGMR 路线: CVODE BDF + SPGMR + 韧性重启
!   29 物种 71 反应 139 命名速率 7 光解
! CLI: mcm_driver [days] [out_pref] [fnox] [fvoc] — 排放缩放 (屋脊实验, 2026-10-25)
! =====================================================================
module mcm_mod
    use oz_mech, only: ns, nr, nq, nj, mre, mpe, update_rates
    implicit none
    integer :: nreact(nr), reac_id(mre, nr)
    integer :: nprod(nr), prod_id(mpe, nr)
    real(8)  :: prod_coef(mpe, nr)
    real(8) :: peff_out(nr), ipr_out(ns, 6), po3c, lo3c
    character(len=8), parameter :: sp_name_arr(29) = &
        [character(len=8) :: "HCHO", "CH3NO3", "CH3OH", "O1D", "O3", "HO2NO2", "NO3", "N2O5", &
         "H2O2", "NO", "NA", "HO2", "NO2", "CH4", "HSO3", "CO", "CL", "O", &
         "HNO3", "SO3", "SO2", "CH3O", "OH", "H2", "HONO", "CH3O2NO2", &
         "CH3OOH", "SA", "CH3O2"]
    character(len=8), parameter :: iprl(6) = ["P", "L", "E", "D", "B", "X"]
    ! MCM 官方光解表槽位 (j, l, m, n) — photolysis-rates_MCM_v331 抄录
    integer, parameter :: NMCMJ = 12
    integer, parameter :: mcm_slot(NMCMJ) = [1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 41, 51]
    real(8), parameter :: mcm_l(NMCMJ) = [6.073d-5, 4.775d-4, 1.041d-5, 1.165d-2, &
        2.485d-2, 1.747d-1, 2.644d-3, 9.312d-7, 4.642d-5, 6.853d-5, 7.649d-6, 1.588d-6]
    real(8), parameter :: mcm_m(NMCMJ) = [1.743d0, 0.298d0, 0.723d0, 0.244d0, &
        0.168d0, 0.155d0, 0.261d0, 1.23d0, 0.762d0, 0.477d0, 0.682d0, 1.154d0]
    real(8), parameter :: mcm_n(NMCMJ) = [0.474d0, 0.08d0, 0.279d0, 0.267d0, &
        0.108d0, 0.125d0, 0.288d0, 0.307d0, 0.353d0, 0.323d0, 0.279d0, 0.318d0]
    real(8), parameter :: PI2 = 3.14159265358979d0
    ! 城市排放 [molec/cm3/s]
    real(8) :: E(ns) = 0.0d0
    ! 干沉降 [1/s]
    real(8) :: vdk(ns) = 0.0d0
    ! 背景松弛
    real(8), parameter :: TAU = 8.64d4
    real(8), parameter :: O3_BG = 6.2e11, CO_BG = 1.5e12, CH4_BG = 4.5e13
    ! 物种索引 (MCM VARIABLE 顺序)
    integer, parameter :: iHCHO=1, iO3=5, iHO2=11, iNO=10, iNO2=13, &
                          iCH4=14, iCO=16, iOH=19, iO=17, iHNO3=22
contains
    subroutine read_stoich(g)
        character(len=*), intent(in) :: g
        integer :: u, k, i
        open(newunit=u, file=g//"/oz.reac", status="old")
        do k = 1, nr
            read(u, *) nreact(k), (reac_id(i, k), i=1, nreact(k))
        end do; close(u)
        open(newunit=u, file=g//"/oz.prod", status="old")
        do k = 1, nr
            read(u, *) nprod(k), (prod_id(i, k), prod_coef(i, k), i=1, nprod(k))
        end do; close(u)
        ! 城市干沉降
        vdk = 0.0d0
        vdk(iO3) = 8.0d-6; vdk(iHNO3) = 2.0d-5; vdk(iHCHO) = 3.0d-6
        ! 城市排放 (日型)
        E(iNO) = 2.0d7      ! NO 20 ppb/h 城市交通
        E(iNO2) = 1.0d6     ! NO2
        E(iCO) = 5.0d6      ! CO
        E(iCH4) = 5.0e3     ! CH4
        E(iHCHO) = 2.0d5    ! HCHO (一次排放)
    end subroutine

    pure function dayf(t)
        real(8), intent(in) :: t
        real(8) :: dayf, tt
        real(8), parameter :: TR = 6.0d0*3600.0d0, TD = 12.0d0*3600.0d0
        tt = mod(t, 86400.0d0)
        dayf = merge(sin(3.14159265358979d0*(tt-TR)/TD)**2, 0.0d0, &
                     tt >= TR .and. tt <= TR+TD)
    end function

    subroutine eval_J_mcm(t, JARR)
        ! SZA(40N, doy172, 时角) → J; 夜间/无表槽位 = 0
        real(8), intent(in) :: t
        real(8), intent(out) :: JARR(nj)
        real(8) :: cosz, decl, lat, ha, hour
        integer :: i
        hour = mod(t, 86400.0d0)/3600.0d0
        decl = 23.45d0*sin(2.0d0*PI2*(172.0d0 - 81.0d0)/365.0d0)*PI2/180.0d0
        lat = 40.0d0*PI2/180.0d0
        ha = 2.0d0*PI2*(hour - 12.0d0)/24.0d0
        cosz = sin(lat)*sin(decl) + cos(lat)*cos(decl)*cos(ha)
        JARR = 0.0d0
        if (cosz <= 1.0d-3) return
        do i = 1, NMCMJ
            JARR(mcm_slot(i)) = mcm_l(i)*cosz**mcm_m(i)*exp(-mcm_n(i)/cosz)
        end do
    end subroutine

    subroutine resid(t, y, ydot)
        real(8), intent(in)  :: t, y(ns)
        real(8), intent(out) :: ydot(ns)
        real(8) :: q(nq), p(nr), JARR(nj), df
        integer :: k, i
        df = dayf(t)
        ! MCM l/m/n 光解 (J = l·cosχ^m·exp(−n/cosχ), lat=40N doy=172) — M4 修复:
        ! 旧版 7 值×df×100 对 nj=51 机制越界 (J<41>/<51> 读垃圾) → 黎明崩溃
        call eval_J_mcm(t, JARR)
        call update_rates(298.0d0, 2.46d19, 5.2d18, 4.45d17, JARR, y, q, p)
        ! 浓度因子 (原始 y — RHS 平滑性优先, 10b#11; kink 曾致 15× 减速)
        do k = 1, nr
            do i = 1, nreact(k)
                p(k) = p(k) * y(reac_id(i, k))
            end do
        end do
        peff_out = p
        ! O3 快收支 (conc csv 两列): P=Σcoefs·peff(产物含 O3), L=Σpeff(反应物含 O3)
        po3c = 0.0d0
        lo3c = 0.0d0
        do k = 1, nr
            do i = 1, nprod(k)
                if (prod_id(i, k) == iO3) po3c = po3c + prod_coef(i, k)*p(k)
            end do
            do i = 1, nreact(k)
                if (reac_id(i, k) == iO3) lo3c = lo3c + p(k)
            end do
        end do
        ! resid: 排放(日型) - 沉降 + 背景 + 化学  (IPR 六账同步分录)
        ydot = 0.0d0
        ipr_out = 0.0d0
        do k = 1, ns
            ydot(k) = E(k) * (0.1d0 + 0.9d0*df) - vdk(k)*y(k)
            ipr_out(k, 3) = E(k) * (0.1d0 + 0.9d0*df)
            ipr_out(k, 4) = -vdk(k)*y(k)
            if (k == iO3) then
                ydot(k) = ydot(k) + (O3_BG - y(k))/TAU
                ipr_out(k, 5) = (O3_BG - y(k))/TAU
            end if
            if (k == iCO) then
                ydot(k) = ydot(k) + (CO_BG - y(k))/TAU
                ipr_out(k, 5) = (CO_BG - y(k))/TAU
            end if
            if (k == iCH4) then
                ydot(k) = ydot(k) + (CH4_BG - y(k))/TAU
                ipr_out(k, 5) = (CH4_BG - y(k))/TAU
            end if
        end do
        do k = 1, nr
            do i = 1, nreact(k)
                ydot(reac_id(i, k)) = ydot(reac_id(i, k)) - p(k)
                ipr_out(reac_id(i, k), 2) = ipr_out(reac_id(i, k), 2) - p(k)
            end do
            do i = 1, nprod(k)
                ydot(prod_id(i, k)) = ydot(prod_id(i, k)) + prod_coef(i, k)*p(k)
                ipr_out(prod_id(i, k), 1) = ipr_out(prod_id(i, k), 1) + prod_coef(i, k)*p(k)
            end do
        end do
    end subroutine
end module

program mcm_run
    use mcm_mod
    implicit none
    real(8) :: Y(ns), T, TOUT, RTOL, ATOL(ns), ROUT(10), RPAR, dum(ns)
    real(8) :: t_end, dt_out, wall0, wall1
    integer(8) :: NEQ8, IOUT(25), IPAR
    integer :: IER, step, nsteps, u1, u2, mxsteps
    character(len=512) :: arg
    character(len=512) :: out_pref = "mcm_city"   ! 长绝对路径需 >256 (旧 128+arg 64 会截断)
    real(8) :: fnox = 1.0d0, fvoc = 1.0d0, days = 5.0d0
    integer :: u3, k, i, islash

    call get_command_argument(1, arg)
    if (len_trim(arg) > 0) read (arg, *) days            ! v2.3.1: 原 arg1/2 被静默忽略
    call get_command_argument(2, arg)
    if (len_trim(arg) > 0) out_pref = trim(arg)
    call get_command_argument(3, arg)
    if (len_trim(arg) > 0) read (arg, *) fnox
    call get_command_argument(4, arg)
    if (len_trim(arg) > 0) read (arg, *) fvoc

    call read_stoich("gen")
    ! 情景缩放 (read_stoich 之后!): NOx 组 (NO/NO2) × fnox, VOC 组 (HCHO) × fvoc
    E(iNO) = E(iNO)*fnox
    E(iNO2) = E(iNO2)*fnox
    E(iHCHO) = E(iHCHO)*fvoc
    Y = 0.0d0
    ! 城市清洁大气初值
    Y(iO3) = 3.0e11; Y(iNO) = 5.0e8; Y(iNO2) = 2.0e8
    Y(iCO) = 2.0e12; Y(iCH4) = 4.5e13; Y(iHCHO) = 2.0e9
    Y(iHO2) = 1.0e8; Y(iOH) = 1.0e6

    RTOL = 1.0d-5
    ATOL = 1.0d4
    ATOL(iO3) = 1.0d9; ATOL(iCH4) = 1.0d11; ATOL(iCO) = 1.0d9
    t_end = days * 86400.0d0
    dt_out = 3600.0d0
    nsteps = int(t_end / dt_out)

    NEQ8 = ns
    call cpu_time(wall0)
    call FNVINITS(1, NEQ8, IER)
    call FCVMALLOC(0.0d0, Y, 2, 2, 2, RTOL, ATOL, IOUT, ROUT, IPAR, RPAR, IER)
    if (IER /= 0) stop 'FCVMALLOC FAIL'
    call FCVDENSE(NEQ8, IER)  ! 29×29 DENSE (NEQ 须 INTEGER*8) (无预条件 SPGMR 黎明崩塌, 10b#12)
    mxsteps = 500000
    call FCVSETIIN('MAX_NSTEPS', mxsteps, IER)
    call FCVSETRIN('MAX_STEP', 60.0d0, IER)

    ! 浓度列单位 = molec/cm3 (2026-10-25 更正: 原 _ppb 表头为存量错标)
    islash = index(trim(out_pref), "/", back=.true.)
    if (islash > 1) call execute_command_line("mkdir -p "//out_pref(1:islash-1))
    open(newunit=u1, file=trim(out_pref)//"_conc.csv", status="replace")
    open(newunit=u2, file=trim(out_pref)//"_irr.csv", status="replace")
    open(newunit=u3, file=trim(out_pref)//"_ipr.csv", status="replace")
    write(u3, '(A)', advance="no") "h"
    do k = 1, ns
        do i = 1, 6
            write(u3, '(",", A, "_", A)', advance="no") trim(sp_name_arr(k)), iprl(i)
        end do
    end do
    write(u3, *)
    write(u1, '(A)') "h,O3,NO,NO2,OH,HO2,HCHO,CO,H2O2,HNO3,NO3,N2O5,HONO,P_O3,L_O3"
    write(u2, '(A)') "h," // repeat("r,", nr-1) // "r"

    do step = 0, nsteps
        if (step > 0) then
            TOUT = dble(step)*dt_out
            call FCVODE(TOUT, T, Y, 1, IER)
            if (IER /= 0) then
                ! 韧性重启
                call FCVREINIT(T, Y, 2, RTOL, ATOL, IER)
                call FCVDENSE(NEQ8, IER)
                call FCVSETRIN('MAX_STEP', 15.0d0, IER)
                call FCVODE(TOUT, T, Y, 1, IER)
                if (IER /= 0) then
                    print *, 'FAIL at h=', step
                    exit
                end if
            end if
        end if
        if (step >= nsteps - 24) then
            call resid(dble(step)*dt_out, Y, dum)
            write(u1, '(I6, 14(",", es13.5))') step, &
                max(Y(iO3),0d0), max(Y(iNO),0d0), max(Y(iNO2),0d0), &
                max(Y(iOH),0d0), max(Y(iHO2),0d0), max(Y(iHCHO),0d0), &
                max(Y(iCO),0d0), max(Y(9),0d0), max(Y(19),0d0), &
                max(Y(7),0d0), max(Y(8),0d0), max(Y(25),0d0), po3c, lo3c
            write(u2, '(I6, 71(",", es13.5))') step, peff_out
            write(u3, '(I6)', advance="no") step
            do k = 1, ns
                do i = 1, 6
                    write(u3, '(",", es13.5)', advance="no") ipr_out(k, i)
                end do
            end do
            write(u3, *)
        end if
    end do
    close(u1); close(u2); close(u3)
    call cpu_time(wall1)
    print '(A, F10.3, A, I0, A, I0)', "MCM-CITY ", wall1-wall0, "s CPU, nst=", IOUT(3), &
        ", netf=", IOUT(5)
end program

subroutine FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
    use mcm_mod, only: resid
    implicit none
    integer(8) :: IPAR(*)
    integer :: IER
    real(8) :: T, Y(*), YDOT(*), RPAR(*)
    call resid(T, Y, YDOT)
    IER = 0
end subroutine
