! =====================================================================
! oz_fac_v3.f90 —— ozone_lab/FAC-V3 工程主力驱动
!   输入: runtime.dat (12 类配置, cfgbuild 编译) — 全功能消费:
!     真光解 (SZA + MCM l/m/n) / 时变气象 (T,P,RH,BLH,RAIN) /
!     时变边界层 (膨胀稀释) / 排放日型 / 干湿沉降 / BC 松弛 / 情景因子
!   求解: CVODE BDF + DENSE (111 物种稠密直接法; 批量并行版同为 DENSE,
!         并行只发生在案例之间, 单箱内部仍单线程)
!   输出: <pref>_conc.csv (111 物种全列) + <pref>_irr.csv (203 反应 IRR)
!         + <pref>_ipr.csv (111 物种 × 6 账 P/L/E/D/B/X 小时窗积分)
! CLI: fac_v3 [scenario_name] [out_prefix] [runtime.dat]
! =====================================================================
module facv3_mod
    use oz_mech, only: ns, nr, nq, nj, mre, mpe, update_rates
    use runtime_dat_mod
    use chem_phys_mod, only: phys_tend, tune_atol_by_name, NPHYS
    implicit none
    integer :: nreact(nr), reac_id(mre, nr)
    integer :: nprod(nr), prod_id(mpe, nr)
    real(8) :: prod_coef(mpe, nr)
    real(8) :: peff_out(nr), ipr_out(ns, 6), ipr_int(ns, 6)
    real(8) :: ipr_a(ns, 6), ipr_w(ns, 6)  ! 梯形 a 端 / 上小时窗累计
    real(8) :: irr_int(nr), irr_w(nr), peff_a(nr)  ! IRR 逐小时积分 (与 IPR 同法)
    character(len=1), parameter :: ipr_lab(6) = ["P", "L", "E", "D", "B", "X"]
contains
    subroutine read_stoich(g)
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
    end subroutine

    subroutine rates_env(t, y, peff, JARR)
        ! 环境驱动的全速率: 气象插值 → M/O2/H2O → update_rates → peff
        real(8), intent(in) :: t, y(ns)
        real(8), intent(out) :: peff(nr), JARR(nj)
        real(8) :: q(nq), p(nr)
        real(8) :: TEMP, PRES, RH, BLH, RAIN, JF, dB, MM, O2, H2O
        integer :: k, i
        call eval_met(t, TEMP, PRES, RH, BLH, RAIN, JF, dB)
        call eval_env(TEMP, PRES, RH, MM, O2, H2O)
        call eval_J(t, JARR)
        call update_rates(TEMP, MM, O2, H2O, JARR, y, q, p)
        peff = p
        do k = 1, nr
            do i = 1, nreact(k)
                peff(k) = peff(k)*y(reac_id(i, k))   ! 原始 y (RHS 平滑性优先)
            end do
        end do
    end subroutine

    subroutine resid_ipr(t, y, ydot)
        ! 六本账: 1机制P 2机制L 3排放 4干+湿沉 5BC松弛+稀释 6BLH膨胀
        real(8), intent(in) :: t, y(ns)
        real(8), intent(out) :: ydot(ns)
        real(8) :: peff(nr), JARR(nj)
        real(8) :: yphys(ns), iphys(ns, NPHYS), metv(7)
        integer :: k, i
        call rates_env(t, y, peff, JARR)
        ydot = 0.0d0
        ipr_out = 0.0d0
        do k = 1, nr
            do i = 1, nreact(k)
                ydot(reac_id(i, k)) = ydot(reac_id(i, k)) - peff(k)
                ipr_out(reac_id(i, k), 2) = ipr_out(reac_id(i, k), 2) - peff(k)
            end do
            do i = 1, nprod(k)
                ydot(prod_id(i, k)) = ydot(prod_id(i, k)) + prod_coef(i, k)*peff(k)
                ipr_out(prod_id(i, k), 1) = ipr_out(prod_id(i, k), 1) + prod_coef(i, k)*peff(k)
            end do
        end do
        ! 物理五项 → chem_phys_mod 唯一实现 (A1; F3 序列 BC 在该模块接线)
        call phys_tend(t, y, emis_scale, yphys, iphys, metv)
        ydot(1:ns) = ydot(1:ns) + yphys(1:ns)
        do k = 1, ns
            ipr_out(k, 3) = iphys(k, 1)
            ipr_out(k, 4) = iphys(k, 2)
            ipr_out(k, 5) = iphys(k, 3)
            ipr_out(k, 6) = iphys(k, 4)
        end do
        peff_out = peff
    end subroutine
end module

program fac_v3
    use facv3_mod
    implicit none
    real(8) :: Y(ns), T, TOUT, RTOL, ATOL(ns), ROUT(10), RPAR, dum(ns)
    real(8) :: t_prev, dt_s
    integer(8) :: NEQ8, IOUT(25), IPAR
    integer :: IER, step, nsteps, nspin, u1, u2, u3, mxsteps, isc, ios, k, j
    character(len=128) :: a1, pref, rtfile
    real(8) :: tstart, tend, o3max
    integer :: iO3v
    call get_command_argument(1, a1)
    call get_command_argument(2, pref)
    call get_command_argument(3, rtfile)
    if (len_trim(rtfile) == 0) rtfile = "runtime_v3.dat"
    if (len_trim(pref) == 0) pref = "output/facv3"
    call ensure_dir(trim(pref))
    call read_runtime(trim(rtfile))
    call read_stoich("gen_v3")
    ! 情景选择 (名或序号)
    isc = 1
    ios = 0
    if (len_trim(a1) > 0) then
        do isc = 1, nscen
            if (trim(scen_name(isc)) == trim(a1)) exit
        end do
        if (isc > nscen) then
            read (a1, *, iostat=ios) isc
            if (ios /= 0 .or. isc > nscen) isc = 1
        end if
    end if
    call scen_build(isc)
    print '(A,A,A,F5.2,A,F5.2,A)', "FAC-V3 情景: ", trim(scen_name(isc)), &
        " (fNOx=", scen_fnox(isc), ", fVOC=", scen_fvoc(isc), ")"
    ! 初值 + 容差
    Y = y0
    RTOL = rtol_cfg
    ATOL = atol_cfg
    call tune_atol_by_name(ATOL)
    NEQ8 = ns
    call cpu_time(tstart)
    call FNVINITS(1, NEQ8, IER)
    call FCVMALLOC(0.0d0, Y, 2, 2, 2, RTOL, ATOL, IOUT, ROUT, IPAR, RPAR, IER)
    call FCVDENSE(NEQ8, IER)
    mxsteps = 200000
    call FCVSETIIN('MAX_NSTEPS', mxsteps, IER)
    call FCVSETRIN('MAX_STEP', 100.0d0, IER)
    iO3v = spidx("O3")
    o3max = 0.0d0
    nsteps = nint(days*86400.0d0/dtout)
    ! spin-up 闸门: step < nspin 只积分不落盘 (spinup_cfg 单位=天, 0=全落盘)
    nspin = nint(spinup_cfg*86400.0d0/dtout)
    if (nspin >= nsteps) then
        print '(A,F4.1,A,F4.1,A)', "FAC-V3 FATAL: spin-up (", spinup_cfg, &
            " d) >= 模拟天数 (", days, " d) → 无任何输出行"
        stop 1
    end if
    open (newunit=u1, file=trim(pref)//"_conc.csv", status="replace")
    open (newunit=u2, file=trim(pref)//"_irr.csv", status="replace")
    open (newunit=u3, file=trim(pref)//"_ipr.csv", status="replace")
    write (u1, '(A)') trim("h"//species_csv_header())
    write (u2, '(A)') trim("h"//irr_csv_header())
    ! IPR 列 = 上一小时窗内的精确时间积分 [molec/cm3] (非快照)
    write (u3, '(A)', advance="no") "h"
    do k = 1, ns
        do j = 1, 6
            write (u3, '(",", A, "_", A)', advance="no") trim(sp_name(k)), &
                ipr_lab(j)
        end do
    end do
    write (u3, *)
    ! ---- ONE_STEP 内步推进: IPR 按内步精确累积 (消除小时快照混叠) ----
    T = 0.0d0
    t_prev = 0.0d0
    ipr_int = 0.0d0
    ipr_w = 0.0d0
    irr_int = 0.0d0
    irr_w = 0.0d0
    call resid_ipr(0.0d0, Y, dum)                ! 梯形 a 端初值 @ (0, Y0)
    ipr_a(1:ns, 1:6) = ipr_out(1:ns, 1:6)
    peff_a(1:nr) = peff_out(1:nr)
    do step = 1, nsteps
        TOUT = dble(step)*dtout
        do while (T < TOUT - 1.0d-9)
            call FCVODE(TOUT, T, Y, 2, IER)          ! ITASK=2: 一个内步即返
            if (IER /= 0) then
                call FCVREINIT(T, Y, 2, RTOL, ATOL, IER)
                call FCVDENSE(NEQ8, IER)
                call FCVSETRIN('MAX_STEP', 10.0d0, IER)
                call FCVODE(TOUT, T, Y, 2, IER)
                if (IER /= 0) then
                    print '(A,I0)', "FAIL @h ", step
                    stop 1
                end if
            end if
            dt_s = min(T, TOUT) - t_prev
            if (dt_s > 0.0d0) then
                call resid_ipr(T, Y, dum)            ! b 端分解 @ (T, Y(T))
                ipr_int(1:ns, 1:6) = ipr_int(1:ns, 1:6) &
                    + 0.5d0*(ipr_a(1:ns, 1:6) + ipr_out(1:ns, 1:6))*dt_s
                ipr_a(1:ns, 1:6) = ipr_out(1:ns, 1:6)
                irr_int(1:nr) = irr_int(1:nr) &
                    + 0.5d0*(peff_a(1:nr) + peff_out(1:nr))*dt_s
                peff_a(1:nr) = peff_out(1:nr)
            end if
            t_prev = T
        end do
        if (T > TOUT + 1.0d-9) then
            ! ONE_STEP 内步过冲 → NORMAL 模式插值回整点 (小时 PA 时间对齐;
            ! CVODE 内部历史连续, 不重复积分)
            call FCVODE(TOUT, T, Y, 1, IER)
            if (IER /= 0) then
                print '(A,I0)', "FAIL(interp) @h ", step
                stop 1
            end if
        end if
        call resid_ipr(TOUT, Y, dum)
        o3max = max(o3max, Y(iO3v))
        if (step >= nspin) then
            write (u1, '(I5)', advance="no") step
            call write_species_row(u1, Y)
            write (u2, '(I5)', advance="no") step
            do k = 1, nr
                write (u2, '(",", es15.5e3)', advance="no") &
                    irr_int(k) - irr_w(k)   ! 本小时窗 IRR 积分
            end do
            write (u2, *)
            write (u3, '(I5)', advance="no") step
            do k = 1, ns
                do j = 1, 6
                    write (u3, '(",", es15.5e3)', advance="no") &
                        ipr_int(k, j) - ipr_w(k, j)   ! 本小时窗积分
                end do
            end do
            write (u3, *)
        end if
        ! 窗口累计推进到本小时末 (写在 diff 之后; 首诊断窗由 step1 起的逐步更新保证)
        ipr_w(1:ns, 1:6) = ipr_int(1:ns, 1:6)
        irr_w(1:nr) = irr_int(1:nr)
    end do
    close (u1); close (u2); close (u3)
    call cpu_time(tend)
    print '(A,F10.2,A,I0,A,F7.2,A)', "FAC-V3 ok (", tend - tstart, " CPU s, nst=", IOUT(3), &
        ", O3max=", o3max/2.46d10, " ppb)"
contains
    subroutine ensure_dir(p)   ! mkdir -p 输出目录 (前缀含路径时)
        character(len=*), intent(in) :: p
        integer :: slash
        slash = index(trim(p), "/", back=.true.)
        if (slash > 1) call execute_command_line("mkdir -p "//p(1:slash-1))
    end subroutine
    integer function spidx(nm)
        character(len=*), intent(in) :: nm
        integer :: i
        spidx = 1
        do i = 1, ns
            if (trim(sp_name(i)) == trim(nm)) then
                spidx = i; return
            end if
        end do
    end function

    function species_csv_header() result(h)
        character(len=1024) :: h
        integer :: i
        h = ""
        do i = 1, ns
            h = trim(h)//","//trim(sp_name(i))
        end do
    end function
    function irr_csv_header() result(h)
        character(len=8*nr) :: h
        integer :: i
        h = ""
        do i = 1, nr
            write (h(len_trim(h)+1:), '(",r",I0)') i
        end do
    end function
    subroutine write_species_row(u, yv)
        integer, intent(in) :: u
        real(8), intent(in) :: yv(ns)
        integer :: i
        do i = 1, ns
            write (u, '(",", es15.5e3)', advance="no") yv(i)
        end do
        write (u, *)
    end subroutine
end program

! IPAR/RPAR 为 fcmix 固定签名 (未用, gfortran 警告可忽略)
subroutine FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
    use facv3_mod, only: resid_ipr
    implicit none
    integer(8) :: IPAR(*)
    integer :: IER
    real(8) :: T, Y(*), YDOT(*), RPAR(*)
    call resid_ipr(T, Y, YDOT)
    IER = 0
end subroutine
