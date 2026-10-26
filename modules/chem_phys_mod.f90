! =====================================================================
! chem_phys_mod.f90 —— 物理五项的唯一实现 (v3 双胞胎共享, A1/A2 统一)
!   排放(日型×情景缩放) / 干+湿沉降 / BC松弛(常量或序列, F3 接线)+稀释 /
!   BLH 膨胀稀释。机制 P/L 由调用方经 oz_mech+计量组装; 本模块只管物理。
!   情景语义 (A2 统一): 调用方传 esc(ns) 排放缩放 + dtv 温度偏移 + jfv 光解因子,
!   fac_v3 来自 runtime.dat 情景表, batch 来自 fnox×fvoc×jfac×dt 网格。
!   IPR 物理四账: 1排放 2干+湿沉 3BC+稀释 4BLH膨胀 (调用方拼 6 账时后移+2)
! =====================================================================
module chem_phys_mod
    use oz_mech, only: ns
    use runtime_dat_mod
    implicit none
    integer, parameter :: NPHYS = 4

contains
    subroutine phys_tend(t, y, esc, yd, ipr, metv)
        ! 注: 物理五项与 T 无关; metv(1) 含 eval_met 内部的情景升温 (dt_scen)
        real(8), intent(in)  :: t, y(ns), esc(ns)
        real(8), intent(out) :: yd(ns), ipr(ns, NPHYS), metv(7)
        real(8) :: TEMP, PRES, RH, BLH, RAIN, JF, dB, kwet, kdil, bcv, shape
        integer :: k
        call eval_met(t, TEMP, PRES, RH, BLH, RAIN, JF, dB)
        metv = (/TEMP, PRES, RH, BLH, RAIN, JF, dB/)
        yd = 0.0d0
        ipr = 0.0d0
        do k = 1, ns
            ! ① 排放 (BOX 日型 / FLAT)
            if (emis0(k) > 0.0d0) then
                if (emis_box(k) == 1) then
                    shape = box_shape(t, emis_amp(k))
                else
                    shape = 1.0d0
                end if
                yd(k) = yd(k) + emis0(k)*esc(k)*shape
                ipr(k, 1) = emis0(k)*esc(k)*shape
            end if
            ! ② 干沉降 (vd/BLH) + 湿沉降 (a·RAIN^b)
            kwet = wet_a(k)*max(RAIN, 0.0d0)**wet_b(k)
            kdil = dry_v(k)*merge(1.0d0, 0.0d0, dry_on(k) == 1)/BLH + kwet
            yd(k) = yd(k) - kdil*y(k)
            ipr(k, 2) = -kdil*y(k)
            ! ③ BC 松弛 (F3: 序列 bc_s 优先, 常量 bc_c 兜底) + 稀释
            if (bc_s(k) > 0 .and. bc_s(k) <= nser) then
                bcv = (ser_interp(bc_s(k), t) - y(k))/tau_mix
            else if (bc_c(k) > 0.0d0) then
                bcv = (bc_c(k) - y(k))/tau_mix
            else
                bcv = 0.0d0
            end if
            bcv = bcv - dil_k*y(k)
            yd(k) = yd(k) + bcv
            ipr(k, 3) = bcv
            ! ④ BLH 膨胀稀释 (仅增长段)
            if (dB > 0.0d0) then
                yd(k) = yd(k) - y(k)*dB/BLH
                ipr(k, 4) = -y(k)*dB/BLH
            end if
        end do
    end subroutine phys_tend

    ! ---- ATOL 按名调优 (单源: fac_v3/batch 双胞胎共用, 防 A1/A2 漂移) ----
    subroutine tune_atol_by_name(a)
        real(8), intent(inout) :: a(ns)
        integer :: i
        do i = 1, ns
            select case (trim(sp_name(i)))
            case ("CH4"); a(i) = 1.0d10
            case ("CO");  a(i) = 1.0d9
            case ("O3");  a(i) = 1.0d8
            case ("HCHO", "PAN", "MPAN", "ISOP", "SO2", "CH3CHO", "MVK", "MACR", &
                  "MGLY", "RONO2", "OLENIT", "TERPNIT", "TERPOOH")
                a(i) = 1.0d7
            case ("DMS", "MSA", "TERP"); a(i) = 1.0d6
            case ("TERPCHO"); a(i) = 1.0d8
            end select
        end do
    end subroutine tune_atol_by_name

    ! ---- dt 覆写版气象 (batch 情景扫描用) ----
    subroutine eval_met_dt(t, dt, TEMP, PRES, RH, BLH, RAIN, JF, dB)
        real(8), intent(in) :: t, dt
        real(8), intent(out) :: TEMP, PRES, RH, BLH, RAIN, JF, dB
        call eval_met(t, TEMP, PRES, RH, BLH, RAIN, JF, dB)
        TEMP = TEMP + dt
    end subroutine eval_met_dt

    ! ---- jfac 覆写版光解 (绝对语义: R11-P0 修复) ----
    ! batch 专用: jf 是「最终光解缩放因子」(已含 config JFAC × CLI 网格因子),
    ! 因此必须先置 jfac_scen=1.0 再调 eval_J, 否则会与 eval_J 内的
    ! `jf = jfac_scen*met_c(7)` 二次相乘 → J ∝ CLI_jfac × JFAC² (旧 bug)
    subroutine eval_J_jf(t, jf_abs, JARR)
        use oz_mech, only: nj
        use runtime_dat_mod, only: eval_J, jfac_scen, met_c, met_s, nser, ser_interp
        real(8), intent(in) :: t, jf_abs
        real(8), intent(out) :: JARR(nj)
        real(8) :: cfg_jf
        ! 取 config 侧的基础 JFAC (常数或序列), 用于把绝对因子归一化掉
        if (met_s(7) > 0 .and. met_s(7) <= nser) then
            cfg_jf = ser_interp(met_s(7), t)
        else
            cfg_jf = met_c(7)
        end if
        jfac_scen = jf_abs/max(cfg_jf, 1.0d-30)   ! 抵消 eval_J 内的 *cfg_jf
        call eval_J(t, JARR)
        jfac_scen = 1.0d0                          ! 复位, 避免污染后续调用
    end subroutine eval_J_jf
end module chem_phys_mod
