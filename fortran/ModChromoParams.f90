module ModChromoParams
  use iso_c_binding,     ONLY: c_float, c_int, c_ptr, c_f_pointer, c_loc
  implicit none
  type, bind(C) :: chromo_params_C
     ! SOLUTION INDEX

     integer(c_int) :: CNI_  ! Ion density
     integer(c_int) :: CNN_  ! Neutral density
     integer(c_int) :: CNV_  ! Ion velocity
     integer(c_int) :: CNU_  ! Neutral velocity
     integer(c_int) :: CEI_  ! Ion energy
     integer(c_int) :: CEN_  ! Neutral energy

     integer(c_int) :: PNI_  ! Ion density (primitive variable)
     integer(c_int) :: PNN_  ! Neutral density (primitive variable)
     integer(c_int) :: PV_   ! Ion velocity (primitive variable)
     integer(c_int) :: PU_   ! Neutral velocity (primitive variable)
     integer(c_int) :: PPI_  ! Ion pressure (primitive variable)
     integer(c_int) :: PPN_  ! Neutral pressure (primitive variable)

     integer(c_int) :: ns          ! Number of grid points in the x-direction
     integer(c_int) :: num_of_eq   ! Number of equations (or number of unknowns)
     integer(c_int) :: num_of_elem ! Total number of elements (ns * num_of_eq)

     ! constants
     real(c_float) :: gammamono  ! Ratio of specific heats (adiabatic index for a monatomic ideal gas)
     real(c_float) :: alpha_p    ! Divergence error propagation parameter (alpha_p)
     real(c_float) :: pi         ! The value of π (pi)

     ! Physical constants
     real(c_float) :: m_i   ! Mass of proton in kilograms (also m0, normalization constant for mass)
     real(c_float) :: m_n   ! Neutral mass
     real(c_float) :: m_e   ! Electron mass
     real(c_float) :: g     ! Gravitational acceleration in normalized units, default is Sun value of 0.27395 km/s^2
     real(c_float) :: mu_0  ! Permeability of free space (magnetic constant)
     real(c_float) :: k_b   ! Boltzmann constant in SI units
     real(c_float) :: e_    ! Elementary charge in Coulombs

     ! Domain normalization constants
     type(c_ptr) :: ds         ! Grid spacing in the x-direction in normalized domain units (array)
     real(c_float) :: CFL      ! Courant-Friedrichs-Lewy condition

     ! stream-aligned quantities (arrays)
     integer(c_int) :: is_allocated = 0
     type(c_ptr) :: dt
     type(c_ptr) :: xn
     type(c_ptr) :: cellX     ! Array of cell positions in X direction
     type(c_ptr) :: cellY     ! Array of cell positions in Y direction
     type(c_ptr) :: cellZ     ! Array of cell positions in Z direction
     type(c_ptr) :: B_imh      ! Magnetic field at i-1/2 in s direction (array)
     type(c_ptr) :: B_iph      ! Magnetic field at i+1/2 in s direction (array)
     type(c_ptr) :: B_i        ! Magnetic field at i in s direction (array)
     type(c_ptr) :: dinvB_ds   ! Inverse B-field gradient in s direction (array)
     type(c_ptr) :: gPotential_imh ! Gravity Potential at face.
     type(c_ptr) :: gPotential_iph ! Gravity Potential at face.
     type(c_ptr) :: outer_bound0
     type(c_ptr) :: outer_bound1
     type(c_ptr) :: inner_bound0
     type(c_ptr) :: inner_bound1
  end type chromo_params_C

  ! Use Fortran style arrays. Params used in WSA code.
  type :: chromo_params_F
     ! Physical parameters
     real(c_float) :: gammamono, alpha_p, pi
     real(c_float) :: m_i, m_n, m_e, g, mu_0, k_b, e_

     ! Cell parameters
     integer(c_int) :: CNI_, CNN_, CNV_, CNU_, CEI_, CEN_
     integer(c_int) :: PNI_, PNN_, PV_, PU_, PPI_, PPN_
     integer(c_int) :: ns, num_of_eq=6, num_of_elem
     real(c_float) :: CFL

     logical :: is_allocated = .false.
     real(c_float), pointer :: Dt_C(:)
     real(c_float), pointer :: LengthSi_G(:)
     real(c_float), pointer :: Coord_DF(:,:)
     real(c_float), pointer :: BSi_F(:)
     real(c_float), pointer :: State_VC(:,:)
     real(c_float), pointer :: gPotential_F(:)
     real(c_float), pointer :: OuterBoundary(:,:)
     real(c_float), pointer :: InnerBoundary(:,:)
  end type chromo_params_F

  interface
     ! extern "C" int Chromo_sub2ind(const int* i, const int* j, const int* k);
     function Chromo_sub2ind(ns, i, j) bind(C, name="Chromo_sub2ind")
       use iso_c_binding,     ONLY: c_int
       integer(c_int), intent(in) :: ns, i, j
       integer(c_int) :: Chromo_sub2ind
     end function Chromo_sub2ind

     ! extern "C" void Chromo_allocate_params_C(struct Chromo_params* pC);
     subroutine Chromo_AllocateParams_C(pC) bind(C, &
          name="Chromo_allocate_params_C")
       import chromo_params_C
       type(chromo_params_C) :: pC
     end subroutine Chromo_AllocateParams_C

     ! extern "C" void Chromo_deallocate_params_C(struct Chromo_params* pC);
     subroutine Chromo_DeallocateParams_C(pC) bind(C, &
          name="Chromo_deallocate_params_C")
       import chromo_params_C
       type(chromo_params_C) :: pC
     end subroutine Chromo_DeallocateParams_C
  end interface

contains
  !============================================================================

  subroutine Chromo_AllocateParams_F(pF)
    type(chromo_params_F), intent(inout) :: pF
    !--------------------------------------------------------------------------
    if(.not. pF%is_allocated) then
       ! All subscripts start from 0
       allocate(pF%Dt_C(0:pF%ns-1))
       allocate(pF%LengthSi_G(0:pF%ns-1))
       allocate(pF%Coord_DF(0:2, 0:pF%ns-1))
       allocate(pF%BSi_F(0:pF%ns))
       allocate(pF%State_VC(0:pF%num_of_eq-1, 0:pF%ns-1))
       allocate(pF%gPotential_F(0:pF%ns))
       allocate(pF%OuterBoundary(0:pF%num_of_eq-1, 0:1))
       allocate(pF%InnerBoundary(0:pF%num_of_eq-1, 0:1))
       pF%is_allocated = .true.
    end if
  end subroutine Chromo_AllocateParams_F
  !============================================================================

  subroutine Chromo_DeallocateParams_F(pF)
    type(chromo_params_F), intent(inout) :: pF
    !--------------------------------------------------------------------------
    if(pF%is_allocated) then
       deallocate(pF%Dt_C)
       deallocate(pF%LengthSi_G)
       deallocate(pF%Coord_DF)
       deallocate(pF%BSi_F)
       deallocate(pF%State_VC)
       deallocate(pF%gPotential_F)
       deallocate(pF%OuterBoundary)
       deallocate(pF%InnerBoundary)
       pF%is_allocated = .false.
    end if
  end subroutine Chromo_DeallocateParams_F
  !============================================================================

  ! This subroutine syncs non-array members of a type(structure).
  ! src=1: pC is the source; src=2: pF is the source.
  subroutine Chromo_SyncParamDim1Var(pC, pF, src)
    type(chromo_params_C) :: pC
    type(chromo_params_F) :: pF
    integer, intent(in) :: src

    !--------------------------------------------------------------------------
    if(src == 1) then
       pF%CNI_ = pC%CNI_
       pF%CNN_ = pC%CNN_
       pF%CNV_ = pC%CNV_
       pF%CNU_ = pC%CNU_
       pF%CEI_ = pC%CEI_
       pF%CEN_ = pC%CEN_
       pF%PNI_ = pC%PNI_
       pF%PNN_ = pC%PNN_
       pF%PV_ = pC%PV_
       pF%PU_ = pC%PU_
       pF%PPI_ = pC%PPI_
       pF%PPN_ = pC%PPN_
       pF%gammamono = pC%gammamono
       pF%alpha_p = pC%alpha_p
       pF%pi = pC%pi
       pF%m_i = pC%m_i
       pF%m_n = pC%m_n
       pF%m_e = pC%m_e
       pF%g = pC%g
       pF%mu_0 = pC%mu_0
       pF%k_b = pC%k_b
       pF%e_ = pC%e_

       pF%ns = pC%ns
       pF%num_of_eq = pC%num_of_eq
       pF%num_of_elem = pC%num_of_elem
       pF%CFL = pC%CFL
    elseif(src == 2) then
       pC%CNI_ = pF%CNI_
       pC%CNN_ = pF%CNN_
       pC%CNV_ = pF%CNV_
       pC%CNU_ = pF%CNU_
       pC%CEI_ = pF%CEI_
       pC%CEN_ = pF%CEN_
       pC%PNI_ = pF%PNI_
       pC%PNN_ = pF%PNN_
       pC%PV_ = pF%PV_
       pC%PU_ = pF%PU_
       pC%PPI_ = pF%PPI_
       pC%PPN_ = pF%PPN_
       pC%gammamono = pF%gammamono
       pC%alpha_p = pF%alpha_p
       pC%pi = pF%pi
       pC%m_i = pF%m_i
       pC%m_n = pF%m_n
       pC%m_e = pF%m_e
       pC%g = pF%g
       pC%mu_0 = pF%mu_0
       pC%k_b = pF%k_b
       pC%e_ = pF%e_

       pC%ns = pF%ns
       pC%num_of_eq = pF%num_of_eq
       pC%num_of_elem = pF%num_of_elem
       pC%CFL = pF%CFL
    end if
  end subroutine Chromo_SyncParamDim1Var
  !============================================================================

  subroutine Chromo_CopyArray_F2C(farray, cptr, size)
    real(c_float), dimension(:), intent(in) :: farray
    type(c_ptr), intent(inout) :: cptr
    integer(c_int), intent(in) :: size

    real(c_float), dimension(:), pointer :: fptr
    !--------------------------------------------------------------------------
    call c_f_pointer(cptr, fptr, [size])
    fptr = farray
  end subroutine Chromo_CopyArray_F2C
  !============================================================================

  ! This subroutine converts xn_I and pC into a single pF.
  subroutine Chromo_ConvertParam_C2F(pC, pF)
    type(chromo_params_C), intent(in) :: pC
    type(chromo_params_F), intent(out) :: pF

    real(c_float), dimension(:), pointer :: dt, xn, ds, cellX, cellY, cellZ
    real(c_float), dimension(:), pointer :: B_imh, B_iph, B_i, dinvB_ds
    real(c_float), dimension(:), pointer :: gPotential_imh, gPotential_iph
    real(c_float), dimension(:), pointer :: outer_bound0, outer_bound1
    real(c_float), dimension(:), pointer :: inner_bound0, inner_bound1

    integer :: i, j

    ! Copy non-array members
    !--------------------------------------------------------------------------
    call Chromo_SyncParamDim1Var(pC, pF, 1)

    ! Transform c_ptr into fortran array
    ! ds is not allocated, but instead points to pC%ds. So no need for deallocation.
    call c_f_pointer(pC%ds, ds, [pC%ns])
    call c_f_pointer(pC%dt, dt, [pC%ns])
    call c_f_pointer(pC%xn, xn, [pC%num_of_elem])
    call c_f_pointer(pC%cellX, cellX, [pC%ns])
    call c_f_pointer(pC%cellY, cellY, [pC%ns])
    call c_f_pointer(pC%cellZ, cellZ, [pC%ns])
    call c_f_pointer(pC%B_imh, B_imh, [pC%ns])
    call c_f_pointer(pC%B_iph, B_iph, [pC%ns])
    call c_f_pointer(pC%B_i, B_i, [pC%ns])
    call c_f_pointer(pC%dinvB_ds, dinvB_ds, [pC%ns])
    call c_f_pointer(pC%gPotential_imh, gPotential_imh, [pC%ns])
    call c_f_pointer(pC%gPotential_iph, gPotential_iph, [pC%ns])
    call c_f_pointer(pC%outer_bound0, outer_bound0, [pC%num_of_eq])
    call c_f_pointer(pC%outer_bound1, outer_bound1, [pC%num_of_eq])
    call c_f_pointer(pC%inner_bound0, inner_bound0, [pC%num_of_eq])
    call c_f_pointer(pC%inner_bound1, inner_bound1, [pC%num_of_eq])

    call Chromo_AllocateParams_F(pF)
    pF%Dt_C = dt
    pF%LengthSi_G = ds
    pF%Coord_DF(0,:) = cellX
    pF%Coord_DF(1,:) = cellY
    pF%Coord_DF(2,:) = cellZ
    pF%BSi_F(0:pF%ns-1) = B_imh
    pF%BSi_F(pF%ns) = B_iph(pF%ns)
    pF%gPotential_F(0:pF%ns-1) = gPotential_imh
    pF%gPotential_F(pF%ns) = gPotential_iph(pF%ns)
    pF%OuterBoundary(:,0) = outer_bound0
    pF%OuterBoundary(:,1) = outer_bound1
    pF%InnerBoundary(:,0) = inner_bound0
    pF%InnerBoundary(:,1) = inner_bound1

    do i = 0, pF%num_of_eq-1
       do j = 0, pF%ns-1
          pF%State_VC(i,j) = xn(Chromo_sub2ind(pF%ns,j,i)+1)
       end do
    end do

  end subroutine Chromo_ConvertParam_C2F
  !============================================================================

  ! This subroutine converts pF into xn_I and pC separately.
  subroutine Chromo_ConvertParam_F2C(pF, pC)
    type(chromo_params_F), intent(in) :: pF
    type(chromo_params_C), intent(out) :: pC

    real(c_float), dimension(:), pointer :: dt, xn, ds, cellX, cellY, cellZ
    real(c_float), dimension(:), pointer :: B_imh, B_iph, B_i, dinvB_ds
    real(c_float), dimension(:), pointer :: gPotential_imh, gPotential_iph
    real(c_float), dimension(:), pointer :: outer_bound0, outer_bound1
    real(c_float), dimension(:), pointer :: inner_bound0, inner_bound1

    integer :: i, j
    ! Copy non-array members
    !--------------------------------------------------------------------------
    call Chromo_SyncParamDim1Var(pC, pF, 2)

    ! Dynamic Arrays
    allocate(ds(pF%ns))
    allocate(dt(pF%ns))
    allocate(xn(pF%num_of_elem))
    allocate(cellX(pF%ns))
    allocate(cellY(pF%ns))
    allocate(cellZ(pF%ns))
    allocate(B_imh(pF%ns))
    allocate(B_iph(pF%ns))
    allocate(B_i(pF%ns))
    allocate(dinvB_ds(pF%ns))
    allocate(gPotential_imh(pF%ns))
    allocate(gPotential_iph(pF%ns))
    allocate(outer_bound0(pF%num_of_eq))
    allocate(outer_bound1(pF%num_of_eq))
    allocate(inner_bound0(pF%num_of_eq))
    allocate(inner_bound1(pF%num_of_eq))

    ! Transform arrays
    dt = pF%Dt_C
    ds = pF%LengthSi_G
    cellX = pF%Coord_DF(0,:)
    cellY = pF%Coord_DF(1,:)
    cellZ = pF%Coord_DF(2,:)
    B_imh = pF%BSi_F(0:pF%ns-1)
    B_iph = pF%BSi_F(1:pF%ns)
    B_i = (B_imh + B_iph)*0.5
    dinvB_ds = (1.0/B_iph - 1.0/B_imh)/ds
    gPotential_imh = pF%gPotential_F(0:pF%ns-1)
    gPotential_iph = pF%gPotential_F(1:pF%ns)
    outer_bound0 = pF%OuterBoundary(:,0)
    outer_bound1 = pF%OuterBoundary(:,1)
    inner_bound0 = pF%innerBoundary(:,0)
    inner_bound1 = pF%innerBoundary(:,1)

    ! Allocate arrays
    call Chromo_AllocateParams_C(pC)

    ! xn indices are (j,i) instead of (i,j)
    do i = 0, pF%num_of_eq-1
       do j = 0, pF%ns-1
          xn(Chromo_sub2ind(pF%ns,j,i)+1) = pF%State_VC(i,j)
       end do
    end do

    ! Copy arrays
    call Chromo_CopyArray_F2C(dt, pC%dt, pC%ns)
    call Chromo_CopyArray_F2C(xn, pC%xn, pC%num_of_elem)
    call Chromo_CopyArray_F2C(ds, pC%ds, pC%ns)
    call Chromo_CopyArray_F2C(cellX, pC%cellX, pC%ns)
    call Chromo_CopyArray_F2C(cellY, pC%cellY, pC%ns)
    call Chromo_CopyArray_F2C(cellZ, pC%cellZ, pC%ns)
    call Chromo_CopyArray_F2C(B_imh, pC%B_imh, pC%ns)
    call Chromo_CopyArray_F2C(B_iph, pC%B_iph, pC%ns)
    call Chromo_CopyArray_F2C(B_i, pC%B_i, pC%ns)
    call Chromo_CopyArray_F2C(dinvB_ds, pC%dinvB_ds, pC%ns)
    call Chromo_CopyArray_F2C(gPotential_imh, pC%gPotential_imh, pC%ns)
    call Chromo_CopyArray_F2C(gPotential_iph, pC%gPotential_iph, pC%ns)
    call Chromo_CopyArray_F2C(outer_bound0, pC%outer_bound0, pC%num_of_eq)
    call Chromo_CopyArray_F2C(outer_bound1, pC%outer_bound1, pC%num_of_eq)
    call Chromo_CopyArray_F2C(inner_bound0, pC%inner_bound0, pC%num_of_eq)
    call Chromo_CopyArray_F2C(inner_bound1, pC%inner_bound1, pC%num_of_eq)

    ! Deallocate memory
    deallocate(ds)
    deallocate(dt)
    deallocate(xn)
    deallocate(cellX)
    deallocate(cellY)
    deallocate(cellZ)
    deallocate(B_imh)
    deallocate(B_iph)
    deallocate(B_i)
    deallocate(dinvB_ds)
    deallocate(gPotential_imh)
    deallocate(gPotential_iph)
    deallocate(outer_bound0)
    deallocate(outer_bound1)
    deallocate(inner_bound0)
    deallocate(inner_bound1)

  end subroutine Chromo_ConvertParam_F2C
  !============================================================================
end module ModChromoParams
!==============================================================================
