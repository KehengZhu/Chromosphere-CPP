!  Copyright (C) 2002 Regents of the University of Michigan,
!  portions used with permission
!  For more information, see http://csem.engin.umich.edu/tools/swmf
module ModChromosphereCPP
  use iso_c_binding,     ONLY: c_float, c_int, c_ptr, c_f_pointer, c_loc
  use ModChromoParams
  implicit none
  interface
     ! extern "C" void Chromo_init_CPP(float* xn, struct Chromo_params* params);
     subroutine Chromo_init_C(pC) bind(C, name="Chromo_init_CPP")
       import chromo_params_C
       type(chromo_params_C), intent(inout) :: pC
     end subroutine Chromo_init_C

     ! extern "C" void Chromo_advance_Euler(float *xn, const float *dt);
     subroutine Chromo_advance_Euler_C(pC) bind(C, &
          name="Chromo_advance_Euler")
       use iso_c_binding,     ONLY: c_float
       import chromo_params_C
       type(chromo_params_C), intent(inout) :: pC
     end subroutine Chromo_advance_Euler_C

     ! extern "C" float Chromo_cal_dt(const float* xn, float* dt);
     subroutine Chromo_cal_dt_C(pC) bind(C, name="Chromo_cal_dt")
       use iso_c_binding,     ONLY: c_float
       import chromo_params_C
       type(chromo_params_C), intent(in) :: pC
     end subroutine Chromo_cal_dt_C

     ! extern "C" void Chromo_prim2cons(float* p_xn);
     subroutine Chromo_prim2cons(p_xn) bind(C, name="Chromo_prim2cons")
       use iso_c_binding,     ONLY: c_ptr
       type(c_ptr), value :: p_xn
     end subroutine Chromo_prim2cons

     ! extern "C" void Chromo_cons2prim(float* c_xn);
     subroutine Chromo_cons2prim(c_xn) bind(C, name="Chromo_cons2prim")
       use iso_c_binding,     ONLY: c_ptr
       type(c_ptr), value :: c_xn
     end subroutine Chromo_cons2prim
  end interface

contains
  !============================================================================
  ! This subroutine is to transform C_ptr to Fortran array
  subroutine Chromo_init(pF)
    type(chromo_params_F), intent(inout) :: pF
    type(chromo_params_C) :: pC
    integer :: i, j

    !--------------------------------------------------------------------------
    call Chromo_AllocateParams_F(pF)

    ! Physical parameters
    pF%gammamono = 5.0 / 3.0                          ! Ratio of specific heats (adiabatic index for a monatomic ideal gas)
    pF%alpha_p = 0.18                                 ! Divergence error propagation parameter (alpha_p)
    pF%pi = 3.14159265358979323846                    ! Value of π (pi) in Fortran
    pF%m_i = 1.6726219e-27                            ! Mass of proton in kilograms (datum::m_p in C++)

    pF%m_n = pF%m_i
    pF%m_e = 9.10938356e-31
    pF%g = 0.27395e3                ! Gravitational acceleration in normalized units, default is Sun value of 0.27395 km/s^2
    pF%mu_0 = 4 * 3.14159265358979323846e-7           ! Permeability of free space (magnetic constant)
    pF%k_b = 1.380649e-23                             ! Boltzmann constant in SI units
    pF%e_ = 1.602176634e-19                           ! Elementary charge in Coulombs
    !  pF%e_ = 1.0                           ! Elementary charge in Coulombs
    !  pF%m_i = 1.0                            ! Mass of proton in kilograms (datum::m_p in C++)
    !  pF%k_b = 1.0                             ! Boltzmann constant in SI units
    !  pF%mu_0 = 1.0           ! Permeability of free space (magnetic constant)
    !  pF%m_n = 1.0

    pF%LengthSi_G = 1.0/pF%ns
    if(pF%CFL < 1e-6) then
       pF%CFL = 0.25                                     ! Courant-Friedrichs-Lewy number
    end if

    ! Solver relevants
    pF%num_of_eq = 6
    pF%num_of_elem = pF%ns * pF%num_of_eq
    pF%BSi_F = 1.0
    call Chromo_ConvertParam_F2C(pF, pC)
    call Chromo_init_C(pC)
    call Chromo_ConvertParam_C2F(pC, pF)

    call Chromo_DeallocateParams_C(pC)
  end subroutine Chromo_init
  !============================================================================

  ! extern "C" void Chromo_advance_Euler(float *xn, const float *dt);
  subroutine Chromo_advance_Euler(pF)
    use ModChromosphereTest,     ONLY: Test_ModelC7_UpdateBC
    type(chromo_params_F), intent(inout) :: pF
    type(chromo_params_C) :: pC

    !--------------------------------------------------------------------------
    ! Update BC
    call Test_ModelC7_UpdateBC(pF)
    call Chromo_ConvertParam_F2C(pF, pC)
    call Chromo_advance_Euler_C(pC)
    ! Update everything in pF
    call Chromo_ConvertParam_C2F(pC, pF)
    call Chromo_DeallocateParams_C(pC)

  end subroutine Chromo_advance_Euler
  !============================================================================

  subroutine Chromo_cal_dt(pF)
    type(chromo_params_F), intent(inout) :: pF

    type(chromo_params_C) :: pC

    !--------------------------------------------------------------------------
    call Chromo_ConvertParam_F2C(pF, pC)
    call Chromo_cal_dt_C(pC)
    ! Update dt in pF
    call Chromo_ConvertParam_C2F(pC, pF)
    call Chromo_DeallocateParams_C(pC)

  end subroutine Chromo_cal_dt
  !============================================================================

end module ModChromosphereCPP
!==============================================================================
