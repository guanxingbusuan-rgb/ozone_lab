! =====================================================================
! cvode_c_mod.f90 —— SUNDIALS CVODE C API 的 ISO_C_BINDING 直连包装
!   动机: fcmix (FCVODE 等) 用 COMMON 块存问题指针 → 单进程仅 1 个
!        并发 CVODE 问题; 直接 C API 每问题独立 mem → OpenMP 线程级并行
!   验证原型: CVodeCreate/Init/SVtolerances/Dense/CVode/Free +
!             N_VNew/GetArray/SetArray/Destroy_Serial (nm -gU 确认)
!   realtype=double, long=64bit (macOS x64)
! =====================================================================
module cvode_c_mod
    use iso_c_binding
    implicit none
    ! CVODE 常数 (cvode.h / cvode_spils.h)
    integer(c_int), parameter :: CV_BDF = 2, CV_NEWTON = 2, CV_NORMAL = 1

    abstract interface
        ! int (*)(realtype t, N_Vector y, N_Vector ydot, void *user_data)
        function cvrhs_proc(t, y, ydot, ud) result(flag) bind(C)
            import c_double, c_ptr, c_int
            real(c_double), value :: t
            type(c_ptr), value :: y, ydot, ud
            integer(c_int) :: flag
        end function
    end interface

    interface
        function c_nvnew(n) result(v) bind(C, name="N_VNew_Serial")
            import c_long, c_ptr
            integer(c_long), value :: n
            type(c_ptr) :: v
        end function
        subroutine c_nvdestroy(v) bind(C, name="N_VDestroy_Serial")
            import c_ptr
            type(c_ptr), value :: v
        end subroutine
        function c_nvgetarr(v) result(p) bind(C, name="N_VGetArrayPointer_Serial")
            import c_ptr
            type(c_ptr), value :: v
            type(c_ptr) :: p
        end function
        subroutine c_nvsetarr(p, v) bind(C, name="N_VSetArrayPointer_Serial")
            import c_ptr
            type(c_ptr), value :: p
            type(c_ptr), value :: v
        end subroutine
        function c_cvcreate(lmm, iter) result(m) bind(C, name="CVodeCreate")
            import c_int, c_ptr
            integer(c_int), value :: lmm, iter
            type(c_ptr) :: m
        end function
        function c_cvinit(m, f, t0, y0) result(flg) bind(C, name="CVodeInit")
            import c_ptr, c_funptr, c_double, c_int
            type(c_ptr), value :: m
            type(c_funptr), value :: f
            real(c_double), value :: t0
            type(c_ptr), value :: y0
            integer(c_int) :: flg
        end function
        function c_cvsvtol(m, rtol, av) result(flg) bind(C, name="CVodeSVtolerances")
            import c_ptr, c_double, c_int
            type(c_ptr), value :: m
            real(c_double), value :: rtol
            type(c_ptr), value :: av
            integer(c_int) :: flg
        end function
        function c_cvdense(m, n) result(flg) bind(C, name="CVDense")
            import c_ptr, c_long, c_int
            type(c_ptr), value :: m
            integer(c_long), value :: n
            integer(c_int) :: flg
        end function
        function c_cvsetmaxsteps(m, mx) result(flg) bind(C, name="CVodeSetMaxNumSteps")
            import c_ptr, c_long, c_int
            type(c_ptr), value :: m
            integer(c_long), value :: mx
            integer(c_int) :: flg
        end function
        function c_cvsetmaxstep(m, h) result(flg) bind(C, name="CVodeSetMaxStep")
            import c_ptr, c_double, c_int
            type(c_ptr), value :: m
            real(c_double), value :: h
            integer(c_int) :: flg
        end function
        function c_cvsetud(m, ud) result(flg) bind(C, name="CVodeSetUserData")
            import c_ptr, c_int
            type(c_ptr), value :: m, ud
            integer(c_int) :: flg
        end function
        function c_cvode(m, tout, yout, tret, itask) result(flg) bind(C, name="CVode")
            import c_ptr, c_double, c_int
            type(c_ptr), value :: m
            real(c_double), value :: tout
            type(c_ptr), value :: yout
            real(c_double) :: tret
            integer(c_int), value :: itask
            integer(c_int) :: flg
        end function
        subroutine c_cvfree(m) bind(C, name="CVodeFree")
            import c_ptr
            type(c_ptr), value :: m
        end subroutine
        function c_cvgetnumsteps(m, p) result(flg) bind(C, name="CVodeGetNumSteps")
            import c_ptr, c_long, c_int
            type(c_ptr), value :: m
            integer(c_long), intent(out) :: p
            integer(c_int) :: flg
        end function
    end interface
end module cvode_c_mod
