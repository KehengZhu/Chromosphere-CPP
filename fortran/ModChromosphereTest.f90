module ModChromosphereTest
  use ModChromoParams
  implicit none
contains
  !============================================================================
  subroutine Test_ShockTube_IC(pF)
    type(chromo_params_F), intent(inout) :: pF

    real :: n_l, n_r, p_l, p_r
    integer :: i

    !--------------------------------------------------------------------------
    n_l = 1.0
    p_l = 1.0
    n_r = 0.125
    p_r = 0.1

    call Chromo_AllocateParams_F(pF)

    do i = 0,pF%ns-1
       if(i > pF%ns/2) then
          pF%State_VC(pF%CNN_, i) = n_r
          pF%State_VC(pF%CEN_, i) = 1.5*p_r
       else
          pF%State_VC(pF%CNN_, i) = n_l
          pF%State_VC(pF%CEN_, i) = 1.5*p_l
       end if
       pF%State_VC(pF%CNI_, i) = 1.0
       pF%State_VC(pF%CEI_, i) = 1.0
    end do
    pF%LengthSi_G = 1.0/pF%ns
  end subroutine Test_ShockTube_IC
  !============================================================================

  subroutine Test_ModelC7_IC(pF)
    type(chromo_params_F), intent(inout) :: pF
    integer :: ns, i
    real, dimension(15, 4) :: arr
    real, dimension(:), allocatable :: height_F, ne_F, nn_F, T_F
    real, dimension(:), allocatable :: length_C, ne_C, nn_C, T_C
    real, dimension(:), allocatable :: length_E, ne_E, nn_E, T_E
    !--------------------------------------------------------------------------
    ns = pF%ns

    ! F for face
    allocate(height_F(-1:ns+3), ne_F(-1:ns+3), nn_F(-1:ns+3), T_F(-1:ns+3))
    ! C for cell
    allocate(length_C(ns), ne_C(ns), nn_C(ns), T_C(ns))
    ! E for extended cell
    allocate(length_E(-1:ns+2), ne_E(-1:ns+2), nn_E(-1:ns+2), T_E(-1:ns+2))

    call Chromo_AllocateParams_F(pF)
    arr = reshape( [ &
         1.003e+03, 1.903e+17, 2.693e+19, 6.225e+03, &
         1.032e+03, 2.021e+17, 2.179e+19, 6.315e+03, &
         1.065e+03, 2.091e+17, 1.722e+19, 6.400e+03, &
         1.101e+03, 2.104e+17, 1.340e+19, 6.474e+03, &
         1.143e+03, 1.995e+17, 1.009e+19, 6.531e+03, &
         1.214e+03, 1.760e+17, 6.322e+18, 6.576e+03, &
         1.299e+03, 1.489e+17, 3.641e+18, 6.598e+03, &
         1.398e+03, 1.438e+17, 1.880e+18, 6.610e+03, &
         1.520e+03, 1.423e+17, 7.956e+17, 6.623e+03, &
         1.617e+03, 1.267e+17, 4.012e+17, 6.633e+03, &
         1.722e+03, 1.027e+17, 1.916e+17, 6.643e+03, &
         1.820e+03, 7.969e+16, 9.843e+16, 6.652e+03, &
         1.894e+03, 6.450e+16, 6.005e+16, 6.660e+03, &
         1.946e+03, 5.526e+16, 4.371e+16, 6.667e+03, &
         1.989e+03, 4.826e+16, 3.447e+16, 6.674e+03 ], &
         shape(arr), order=[2,1])

    do i=-1, ns+3
       height_F(i) = arr(1,1) + 1.0*(arr(15,1)-arr(1,1))/(size(height_F)-1.0)*(i+1.0)
    end do

    call interp1(arr(:,1), arr(:,2), height_F, ne_F)
    call interp1(arr(:,1), arr(:,3), height_F, nn_F)
    call interp1(arr(:,1), arr(:,4), height_F, T_F)

    length_E = height_F(0:ns+3) - height_F(-1:ns+2)
    length_E = length_E*1000.0
    ne_E = 0.5*(ne_F(0:ns+3) + ne_F(-1:ns+2))
    nn_E = 0.5*(nn_F(0:ns+3) + nn_F(-1:ns+2))
    T_E = 0.5*(T_F(0:ns+3) + T_F(-1:ns+2))

    length_C = length_E(1:ns)
    ne_C = ne_E(1:ns)
    nn_C = nn_E(1:ns)
    T_C = T_E(1:ns)

    pF%gPotential_F = 0
    pF%LengthSi_G = length_C
    pF%BSi_F = 1.0

    pF%State_VC(pF%CNI_, :) = ne_C*pF%m_i
    pF%State_VC(pF%CNN_, :) = nn_C*pF%m_n
    pF%State_VC(pF%CNU_, :) = 0
    pF%State_VC(pF%CNV_, :) = 0
    pF%State_VC(pF%CEI_, :) = 3.0/2.0*pF%k_b*ne_C*2.0*T_C + ne_C*pF%m_i*0.5*(pF%gPotential_F(1:ns)+pF%gPotential_F(0:ns-1))
    pF%State_VC(pF%CEN_, :) = 3.0/2.0*pF%k_b*nn_C*T_C + nn_C*pF%m_n*0.5*(pF%gPotential_F(1:ns)+pF%gPotential_F(0:ns-1))

    ! Outer BC
    do i=0, 1
      pF%OuterBoundary(pF%CNI_, i) = ne_E(ns+1+i)*pF%m_i
      pF%OuterBoundary(pF%CNN_, i) = nn_E(ns+1+i)*pF%m_n
      pF%OuterBoundary(pF%CNV_, i) = 0
      pF%OuterBoundary(pF%CNU_, i) = 0
      pF%OuterBoundary(pF%CEI_, i) = 3.0/2.0*pF%k_b*ne_E(ns+1+i)*(4.0*T_E(ns+1+i)) + ne_E(ns+1+i)*pF%m_i*pF%gPotential_F(ns)
      pF%OuterBoundary(pF%CEN_, i) = 3.0/2.0*pF%k_b*nn_E(ns+1+i)*(2.0*T_E(ns+1+i)) + nn_E(ns+1+i)*pF%m_i*pF%gPotential_F(ns)
    end do

    ! Inner BC
    do i=-1, 0 ! InnerBoundary(1) is more inner than InnerBoundary(0)
      pF%InnerBoundary(pF%CNI_, -i) = ne_E(i)*pF%m_i
      pF%InnerBoundary(pF%CNN_, -i) = nn_E(i)*pF%m_n
      pF%InnerBoundary(pF%CNV_, -i) = 0
      pF%InnerBoundary(pF%CNU_, -i) = 0
      pF%InnerBoundary(pF%CEI_, -i) = 3.0/2.0*pF%k_b*ne_E(i)*2.0*T_E(i) + ne_E(i)*pF%m_i*pF%gPotential_F(1)
      pF%InnerBoundary(pF%CEN_, -i) = 3.0/2.0*pF%k_b*nn_E(i)*T_E(i) + nn_E(i)*pF%m_n*pF%gPotential_F(1)
    end do

    deallocate(height_F, ne_F, nn_F, T_F)
    deallocate(length_C, ne_C, nn_C, T_C)
    deallocate(length_E, ne_E, nn_E, T_E)
  end subroutine Test_ModelC7_IC
  !============================================================================

  subroutine Test_ModelC7_UpdateBC(pF)
    type(chromo_params_F) :: pF
    integer :: ns
    real, dimension(:), allocatable :: ni, nn, Ti, Tn, ei, en, uu, vv
    real, dimension(0:1) :: Ti_B, Tn_B, uu_B, vv_B, ni_B, nn_B

    !--------------------------------------------------------------------------
    ns = pF%ns
    allocate(ni(ns), nn(ns), Ti(ns), Tn(ns), ei(ns), en(ns), uu(ns), vv(ns))

    ni = pF%State_VC(pF%CNI_, :)/pF%m_i
    nn = pF%State_VC(pF%CNN_, :)/pF%m_n
    ei = pF%State_VC(pF%CEI_, :)
    en = pF%State_VC(pF%CEN_, :)
    vv = pF%State_VC(pF%CNV_, :)/(pF%m_i*ni)
    uu = pF%State_VC(pF%CNU_, :)/(pF%m_n*nn)
    Ti = (2.0/3.0*ei - 1.0/3.0*pF%m_i*ni*vv*vv - 2.0/3.0*pF%m_i*ni*0.5*(pF%gPotential_F(0:ns-1)+pF%gPotential_F(1:ns)))/(2.0*ni*pF%k_b)
    Tn = (2.0/3.0*en - 1.0/3.0*pF%m_n*nn*uu*uu - 2.0/3.0*pF%m_n*nn*0.5*(pF%gPotential_F(0:ns-1)+pF%gPotential_F(1:ns)))/(nn*pF%k_b)

    ni_B = pF%OuterBoundary(pF%CNI_, :)/pF%m_i
    nn_B = pF%OuterBoundary(pF%CNN_, :)/pF%m_n
    vv_B = pF%OuterBoundary(pF%CNV_, :)/(pF%m_i*ni_B)
    uu_B = pF%OuterBoundary(pF%CNU_, :)/(pF%m_n*nn_B)
    Ti_B = (2.0/3.0*pF%OuterBoundary(pF%CEI_, :) - 1.0/3.0*pF%m_i*ni_B*vv_B*vv_B - 2.0/3.0*pF%m_i*ni_B*pF%gPotential_F(ns))/(2.0*ni_B*pF%k_b)
    Tn_B = (2.0/3.0*pF%OuterBoundary(pF%CEN_, :) - 1.0/3.0*pF%m_n*nn_B*uu_B*uu_B - 2.0/3.0*pF%m_n*nn_B*pF%gPotential_F(ns))/(nn_B*pF%k_b)

    ! BC for ni, nn, Ti, Tn do not change. uu and vv are halved.
    vv_B(0) = 0.5*vv(ns)
    uu_B(0) = 0.5*uu(ns)
    vv_B(1) = 0.5*vv_B(0)
    uu_B(1) = 0.5*uu_B(0)
    pF%OuterBoundary(pF%CNV_, :) = pF%m_i*ni_B*vv_B
    pF%OuterBoundary(pF%CNU_, :) = pF%m_n*nn_B*uu_B
    pF%OuterBoundary(pF%CEI_, :) = 3.0/2.0*pF%k_b*ni_B*2.0*Ti_B + 0.5*pF%m_i*ni_B*vv_B*vv_B + pF%m_i*ni_B*pF%gPotential_F(ns)
    pF%OuterBoundary(pF%CEN_, :) = 3.0/2.0*pF%k_b*nn_B*Tn_B + 0.5*pF%m_n*nn_B*uu_B*uu_B + pF%m_n*nn_B*pF%gPotential_F(ns)

    deallocate(ni, nn, Ti, Tn, ei, en, uu, vv)
  end subroutine Test_ModelC7_UpdateBC
  !============================================================================
  subroutine interp1(xData, yData, xVal, yVal)
   implicit none
   ! Inputs
   real, intent(in) :: xData(:)       ! x-values of the data to be interpolated
   real, intent(in) :: yData(:)       ! y-values of the data to be interpolated
   real, intent(in) :: xVal(:)        ! x-values where interpolation should be performed
   ! Output
   real, intent(out) :: yVal(size(xVal)) ! interpolated y-values at xVal points
   
   ! Local variables
   integer :: n, i, j
   real, allocatable :: h(:), a(:), b(:), c(:), d(:), alpha(:), l(:), mu(:), z(:), c2(:)
   
   n = size(xData)
   
   if (n < 2) then
       print*, 'Error: At least two data points are required for spline interpolation.'
       stop
   endif

   ! Allocate arrays
   allocate(h(n-1), a(n), b(n-1), c(n), d(n-1), alpha(n), l(n), mu(n), z(n), c2(n-1))

   ! Step 1: Set up the equations
   do i = 1, n-1
       h(i) = xData(i+1) - xData(i)
   end do

   do i = 2, n-1
       alpha(i) = (3.0 / h(i)) * (yData(i+1) - yData(i)) - (3.0 / h(i-1)) * (yData(i) - yData(i-1))
   end do

   ! Step 2: Solve the tridiagonal system
   l(1) = 1.0
   mu(1) = 0.0
   z(1) = 0.0

   do i = 2, n-1
       l(i) = 2.0 * (xData(i+1) - xData(i-1)) - h(i-1) * mu(i-1)
       mu(i) = h(i) / l(i)
       z(i) = (alpha(i) - h(i-1) * z(i-1)) / l(i)
   end do

   l(n) = 1.0
   z(n) = 0.0
   c(n) = 0.0

   do j = n-1, 1, -1
       c(j) = z(j) - mu(j) * c(j+1)
       b(j) = (yData(j+1) - yData(j)) / h(j) - h(j) * (c(j+1) + 2.0 * c(j)) / 3.0
       d(j) = (c(j+1) - c(j)) / (3.0 * h(j))
       a(j) = yData(j)
   end do

   ! Step 3: Perform the interpolation
   do i = 1, size(xVal)
       ! Find the right interval
       j = 1
       do while (j < n .and. xVal(i) > xData(j+1))
           j = j + 1
       end do
       
       yVal(i) = a(j) + b(j) * (xVal(i) - xData(j)) + c(j) * (xVal(i) - xData(j))**2 + d(j) * (xVal(i) - xData(j))**3
   end do

   ! Deallocate arrays
   deallocate(h, a, b, c, d, alpha, l, mu, z, c2)

end subroutine interp1


end module ModChromosphereTest
!==============================================================================
