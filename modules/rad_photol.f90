! =====================================================================
! rad_photol.f90 —— B2: 在线光解 (Day 26 辐射柱自包含, 无交叉引用)
!   计算 J_NO2/J_O3D/J_O3P/J_HCHO/J_H2O2 基于:
!     太阳天顶角(纬度+日序) + O3柱(DU) + 云光学厚度
!   替代 MCM l/m/n 参数化表 → 辐射-化学在线耦合
! CLI: rad_photol lat doy cloud_tau → 输出 7 个 J 值逐时
! =====================================================================
module rad_photol_mod
    implicit none
    real(8), parameter :: PI = 3.14159265358979d0
    real(8), parameter :: S0 = 1361.0d0    ! 太阳常数 W/m²
contains
    function solar_zenith(lat_deg, doy, hour)
        ! 太阳天顶角余弦 (简化: 无均时差)
        real(8) :: solar_zenith, lat_deg, doy, hour
        real(8) :: decl, lat, ha
        decl = 23.45d0 * sin(2.0d0*PI*(doy - 81.0d0)/365.0d0) * PI/180.0d0
        lat = lat_deg * PI/180.0d0
        ha = 2.0d0*PI*(hour - 12.0d0)/24.0d0
        solar_zenith = sin(lat)*sin(decl) + cos(lat)*cos(decl)*cos(ha)
    end function

    function beer_lambert(I0, tau)
        real(8) :: beer_lambert, I0, tau
        beer_lambert = I0 * exp(-tau)
    end function

    subroutine calc_J(lat_deg, doy, hour, o3_du, cloud_tau, JARR)
        ! 在线计算 7 个光解速率
        real(8), intent(in)  :: lat_deg, doy, hour, o3_du, cloud_tau
        real(8), intent(out) :: JARR(7)
        real(8) :: mu, tau_uv, tau_vis, f_dir, f_diff

        mu = solar_zenith(lat_deg, doy, hour)
        if (mu <= 0.0d0) then
            JARR = 0.0d0
            return
        end if

        ! O3 UV 吸收光学厚度 (Hartley 带, 简化)
        ! 300 DU ≈ 8.03e18 molec/cm2; σ(Hartley 300nm) ≈ 1e-19 cm2
        !   → tau(mu=1) ≈ 0.8; 取 0.01/DU 使 300 DU → tau=3 (含斜程增强的保守值)
        ! (旧实现把 o3_du 直接当 tau_uv≈300 → O1D 光解被压到 1e-65, 严重失真)
        tau_uv = 0.01d0 * o3_du * (1.0d0/mu)
        ! Rayleigh + 云 (可见-近红外)
        tau_vis = 0.15d0/mu + cloud_tau
        ! 云透射因子
        f_dir = exp(-cloud_tau/mu)
        f_diff = 0.3d0*(1.0d0 - f_dir)  ! 简化: 30% 漫射

        ! J 值 (Beer-Lambert + 光子通量 → 速率)
        ! 基线取自 MCM 正午典型值, 按 (太阳角×O3×云) 缩放
        JARR(1) = 8.0d-3 * mu * f_dir                    ! J_NO2 (可见, 不受 O3)
        JARR(2) = 2.5d-5 * mu * f_dir * exp(-tau_uv*0.5d0) ! J_O3D (UV, 强 O3 吸收)
        JARR(3) = 2.0d-4 * mu * f_dir * exp(-tau_uv*0.3d0) ! J_O3P
        JARR(4) = 2.8d-4 * mu * (f_dir + f_diff)           ! J_HCHO(radical)
        JARR(5) = 5.0d-4 * mu * (f_dir + f_diff)           ! J_H2O2
        JARR(6) = 1.5d-6 * mu * f_dir * exp(-tau_uv*0.8d0) ! J_HNO3 (深 UV)
        JARR(7) = 3.0d-4 * mu * (f_dir + f_diff)           ! J_MOOH
    end subroutine
end module

program rad_photol_demo
    use rad_photol_mod
    implicit none
    real(8) :: JARR(7), lat, doy, hour, o3, cld
    integer :: h
    lat = 45.0d0; doy = 172.0d0; o3 = 300.0d0; cld = 0.0d0
    print '(A)', "B2 在线光解演示 (45°N, 夏至, 300 DU O3, 无云)"
    print '(A)', "  hour   J_NO2     J_O3D     J_HCHO    J_H2O2"
    do h = 0, 24, 3
        hour = dble(h)
        call calc_J(lat, doy, hour, o3, cld, JARR)
        print '(I6, 4(1x,ES12.3))', h, JARR(1), JARR(2), JARR(4), JARR(5)
    end do
    ! 云敏感性
    print '(A)', "\n云敏感性 (正午):"
    do cld = 0.0d0, 10.0d0, 2.5d0
        call calc_J(lat, doy, 12.0d0, o3, cld, JARR)
        print '("  tau=",F5.1,"  J_NO2=",ES12.3,"  J_O3D=",ES12.3)', cld, JARR(1), JARR(2)
    end do
end program
