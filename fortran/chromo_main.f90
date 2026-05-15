program Chromo_main
  use ModChromosphereCPP
  use ModChromosphereTest
  implicit none

  integer(c_int), parameter :: ns = 100, n_eq = 6
  real(c_float) :: time, V, total_time
  type(chromo_params_F) :: tParams
  type(chromo_params_C) :: pC
  real(c_float), pointer :: xn_I(:), Cs_I(:)
  integer :: i, j, test_id
  character(len=100) :: filename

  ! Initialize the chromosphere cpp code
  !----------------------------------------------------------------------------

  ! Initialize tParams
  tParams%ns = ns
  tParams%CFL = 0.25
  call Chromo_init(tParams)
  ! ShockTube_IC
   !  call Test_ShockTube_IC(tParams)
  call Test_ModelC7_IC(tParams)
  allocate(Cs_I(ns))
  Cs_I = 2e4

  ! Start main loop
  time = 0.0
  open(20, file='output.log', status='unknown')
  total_time = 10*sum(tParams%LengthSi_G)/2e4
  write(20,*) "Total time = ", total_time
  write(*, *) "Total time = ", total_time

  ! Main loop
  i = 0
  do while (time < total_time .and. i < 10000)
     call Chromo_cal_dt(tParams)
     Cs_I = tParams%CFL*tParams%LengthSi_G/tParams%Dt_C
   !   print*, 'dt = ', sum(tParams%Dt_C)/size(tParams%Dt_C)
   !   print*, "cs = ", Cs_I
     if(mod(i, 100) == 0) then
        write(*, *) "step = ", i, "dt = ", sum(tParams%Dt_C)/size(tParams%Dt_C), "time = ", time
     end if
     call Chromo_advance_Euler(tParams)
     time = time + sum(tParams%Dt_C)/size(tParams%Dt_C)
     i = i + 1
  end do

  ! Output
  call output_data(tParams, 'output.txt')
  close(20)

  !   test_id = 0
  !   if(test_id == 0) then
  !      ! Shock tube test
  !      call Chromo_init(xn_I, ns, tParams)
  !      time = 0.0

  !      open(20, file='output.log', status='unknown')
  !      total_time = 0.15
  !      write(20,*) "Total time = ", total_time
  !      write(*, *) "Total time = ", total_time

  !      ! Main loop
  !      i = 0
  !      do while (time < total_time .and. i < 10000)
  !         call Chromo_cal_dt(xn_I, dt)
  !         if(mod(i, 100) == 0) then
  !            write(*, *) "step = ", i, "dt = ", sum(dt)/size(dt), "time = ", time
  !         end if
  !         call Chromo_advance_Euler(xn_I, dt)
  !         time = time + sum(dt)/size(dt)
  !         i = i + 1
  !      end do

  !      ! save to file (temporary text, will use hdf5 later)
  !      open(10, file='output.txt', status='unknown')
  !      call Chromo_cons2prim(xn_I)
  !      do i = 1, ns*n_eq
  !         write(10,*) xn_I(i)
  !      end do
  !      close(10)
  !      write(20, *) "END!"
  !      close(20)
  !   elseif(test_id == 1) then
  !      call Chromo_init(xn_I, ns, tParams)
  !      time = 0

  !      total_time = 0.15
  !      write(*, *) "Total time = ", total_time

  !      ! Main loop
  !      i = 0
  !      print *, "Start main loop"
  !      do while (time < total_time .and. i < 10000)
  !         call Chromo_cal_dt(xn_I, dt)
  !         if(mod(i, 100) == 0) then
  !            write(*, *) "step = ", i, "dt = ", sum(dt)/size(dt), "time = ", time
  !         end if
  !         call Chromo_advance_Euler(xn_I, dt)
  !         time = time + sum(dt)/size(dt)
  !         i = i + 1
  !      end do
  !      print *, "End main loop"

  !      ! save to file (temporary text, will use hdf5 later)
  !      write(filename, '(A,I0,A)') 'output_', ns, '.txt'
  !      open(10, file=filename, status='unknown')
  !      do i = 1, ns*n_eq
  !         write(10,*) xn_I(i)
  !      end do
  !      close(10)

  !      write(*, *) "END!"
  !   endif
contains
   subroutine output_data(pF, file_path)
      type(chromo_params_F) :: pF
      character(len=*) :: file_path
      open(10, file=file_path, status='unknown')
      do i=0, ns-1
         write(10, '(E15.3E3)', advance='no') sum(tParams%LengthSi_G(0:i)) + 1.003e3
      end do
      write(10, *) ''
      do i=0, ns-1
         write(10, *) tParams%State_VC(:, i)
      end do
      close(10)
   end subroutine
end program Chromo_main
!==============================================================================
