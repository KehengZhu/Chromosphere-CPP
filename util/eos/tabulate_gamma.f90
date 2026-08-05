!  Tabulate the classical, pure-hydrogen CRASH equation of state used to
!  validate Chromosphere2026's analytic Saha closure and supply Gamma_1.
!
!  For each material we sweep (temperature T, heavy-particle number density Na)
!  and, at every grid point, solve the Saha ionization equilibrium and read off
!  the CRASH thermodynamics:
!
!     Z      = <Z>            average ionization degree            (CRASH zAv)
!     P      = pressure                                    [Pa]
!     Edens  = internal energy density                     [J/m^3]
!     Gamma  = 1 + P/Edens    ENERGY gamma, i.e. e = P/(Gamma-1)  (get_gamma GammaOut)
!     GammaS = (dlnP/dlnrho)_S ADIABATIC / sound-speed gamma      (get_gamma GammaSOut)
!     Gammae, GammaSe          the electron-only analogues
!     Cv     = heat capacity at constant volume            [per atom, units k_B]
!
!  The default production mode keeps excitation, Fermi-gas, and Coulomb
!  corrections OFF.  The optional "excitation" diagnostic mode enables bound
!  excitation while keeping Fermi and Coulomb corrections off, and writes only
!  under outputs/eos_gamma/ so it cannot replace the signed-off runtime table.
!  Ground-state statistical weights are ON in both modes.
!
!  Not included in gamma: radiative cooling, radiation energy, conduction,
!  photoionization, or finite-rate/non-equilibrium ionization and recombination.
!  Do not combine this equilibrium-ionization energy closure unchanged with a
!  separate source stage that also pays/returns the ionization potential.
!
!  Two files are written for the one pure-H material:
!    gamma_hydrogen.dat  full 11-column validation table
!    data/eos/gamma1_hydrogen_v1.dat runtime table:
!       log10(T), log10(n_H), Gamma_1

program tabulate_gamma

  use CRASH_ModStatSum
  use CRASH_ModPartition,  ONLY: UseCoulombCorrection, zAv, ToleranceZ, &
       StatSumToleranceLog
  use CRASH_ModExcitation, ONLY: UseGroundStatWeight, LogGi_II
  use CRASH_ModExcitationData, ONLY: UseExcitation
  use CRASH_ModFermiGas,   ONLY: UseFermiGas, LogGeMinFermi, LogGeMinBoltzmann
  use CRASH_ModIonization, ONLY: init_ioniz_potential, put_ioniz_potential
  use CRASH_ModAtomicMass, ONLY: cAtomicMass_I
  use ModConst,            ONLY: cEV, cAtomicMass

  implicit none

  ! ---- grid ---------------------------------------------------------------
  integer :: nT, nNa                  ! temperature/density points (log spaced)
  ! CRASH hard-zeros ionization below 0.02 chi_H/k_B (~3157 K), so the lower
  ! edge is placed just above that implementation cutoff.  It remains below
  ! the coolest expected C7/runtime state (~4400 K); the upper edge covers
  ! flare evaporation with a full decade of margin above 10^7 K.
  real,    parameter :: TMin = 3.2e3, TMax = 1.0e8           ! [K]
  real,    parameter :: LogNaMin = 12.0, LogNaMax = 26.0     ! log10(n_H [m^-3])

  ! ---- material description ----------------------------------------------
  character(len=256) :: FileValidation, FileRuntime
  character(len=32)  :: Mode
  logical :: IsExcitation

  integer :: iT, iNa, iError, iUnit, iRuntime
  real    :: T_K, TeEV, NaSi, LogNa, Rho, Amean
  real    :: P, Edens, Gamma, GammaS, Gammae, GammaSe, Cv, ZbarOut
  real    :: dLogT
  real, parameter :: ChiH_J = 2.179872361e-18
  real, parameter :: BoltzmannSI = 1.380649e-23

  !-------------------------------------------------------------------------
  ! Classical pure-H Saha model.  Keep these assignments explicit so a change
  ! in CRASH library defaults cannot silently change the generated table.
  UseExcitation         = .false.
  UseFermiGas           = .false.
  UseCoulombCorrection  = .false.
  UseGroundStatWeight   = .true.
  LogGeMinFermi         = -4.0
  ! With UseFermiGas=.false. CRASH otherwise still clamps LogGe at its default
  ! value 4.  A very low floor keeps the requested model genuinely Boltzmann.
  LogGeMinBoltzmann     = -700.0
  ToleranceZ            = 1.0e-12
  StatSumToleranceLog   = 700.0
  UsePreviousTe        = .true.

  ! Default generation remains the tracked 501x57 production table.  The
  ! optional refined mode doubles both interval counts, preserving every
  ! production node and adding direct-CRASH midpoint calls for Stage-9 tests.
  nT = 501
  nNa = 57
  IsExcitation = .false.
  FileValidation = 'outputs/eos_gamma/gamma_hydrogen.dat'
  FileRuntime = 'data/eos/gamma1_hydrogen_v1.dat'
  Mode = ''
  call get_command_argument(1, Mode)
  if(trim(Mode) == 'refined')then
     nT = 1001
     nNa = 113
     FileValidation = 'outputs/eos_gamma/gamma_hydrogen_refined.dat'
     FileRuntime = 'outputs/eos_gamma/gamma1_hydrogen_refined.dat'
  else if(trim(Mode) == 'excitation')then
     IsExcitation = .true.
     UseExcitation = .true.
     FileValidation = 'outputs/eos_gamma/gamma_hydrogen_excitation.dat'
     FileRuntime = 'outputs/eos_gamma/gamma1_hydrogen_excitation.dat'
  else if(len_trim(Mode) /= 0)then
     write(*,'(a)') 'usage: tabulate_gamma.exe [refined|excitation]'
     error stop 2
  end if

  dLogT = log(TMax/TMin) / real(nT-1)

  ! Match Chromosphere2026's existing hydrogen energy zero exactly.  Set the
  ! database value before set_element copies it into the active material.
  call init_ioniz_potential
  call put_ioniz_potential(1, 1, ChiH_J/cEV)
  call set_element(1)
  ! With excitation disabled, apply the neutral-H ground degeneracy explicitly.
  ! In excitation mode CRASH reads the complete bound-level degeneracies.
  if(.not.IsExcitation) LogGi_II(0,1) = log(2.0)
  Amean = cAtomicMass_I(1)

  iUnit = 20
  iRuntime = 21
  open(iUnit, file=trim(FileValidation), status='replace')
  open(iRuntime, file=trim(FileRuntime), status='replace')
     if(IsExcitation)then
        write(iUnit,'(a)') '# CRASH pure-H Saha-Boltzmann EOS -- excitation diagnostic table'
        write(iUnit,'(a)') '# excitation=ON Fermi=OFF Coulomb=OFF ground-stat-weight=ON'
     else
        write(iUnit,'(a)') '# CRASH classical pure-H Saha EOS -- validation table'
        write(iUnit,'(a)') '# excitation=OFF Fermi=OFF Coulomb=OFF ground-stat-weight=ON'
     end if
     write(iUnit,'(a)') '# Gamma  = 1 + P/Edens        (energy gamma, closes e = P/(Gamma-1))'
     write(iUnit,'(a)') '# GammaS = (dlnP/dlnrho)_S     (adiabatic / sound-speed gamma)'
     write(iUnit,'(a)') '# columns:'
     write(iUnit,'(a)') '#   1 T[K]  2 Na[m^-3]  3 Rho[kg/m^3]  4 Zbar  5 P[Pa]  '// &
          '6 Edens[J/m^3]  7 Gamma  8 GammaS  9 Gammae  10 GammaSe  11 Cv[k_B/atom]'
     write(iRuntime,'(a)') '# format_version=1'
     if(IsExcitation)then
        write(iRuntime,'(a)') '# table_id=chromosphere2026_gamma1_hydrogen_excitation'
     else
        write(iRuntime,'(a)') '# table_id=chromosphere2026_gamma1_hydrogen_v1'
     end if
     write(iRuntime,'(a)') '# material=pure_H'
     if(IsExcitation)then
        write(iRuntime,'(a)') '# excitation=1'
     else
        write(iRuntime,'(a)') '# excitation=0'
     end if
     write(iRuntime,'(a)') '# fermi_gas=0'
     write(iRuntime,'(a)') '# coulomb_correction=0'
     write(iRuntime,'(a)') '# ground_stat_weight=1'
     if(IsExcitation)then
        write(iRuntime,'(a)') '# saha_prefactor=bound_partition_function'
     else
        write(iRuntime,'(a)') '# saha_prefactor=coefficient_one'
     end if
     write(iRuntime,'(a)') '# axis_1=log10_T_K'
     write(iRuntime,'(a)') '# axis_2=log10_nH_m-3'
     write(iRuntime,'(a)') '# value=Gamma1'
     write(iRuntime,'(a,i0)') '# nT=', nT
     write(iRuntime,'(a,i0)') '# nN=', nNa
     write(iRuntime,'(a)') '# k_b_J_K=1.380649e-23'
     write(iRuntime,'(a)') '# m_e_kg=9.1093837015e-31'
     write(iRuntime,'(a)') '# m_H_kg=1.6726219e-27'
     write(iRuntime,'(a)') '# h_J_s=6.62607015e-34'
     write(iRuntime,'(a)') '# chi_H_J=2.179872361e-18'
     write(iRuntime,'(a)') '# generator=util/eos/tabulate_gamma.f90'
     write(iRuntime,'(a)') '# generation_date=2026-07-19'
     write(iRuntime,'(a)') '# source_revision=crash_gamma_eos_stage0_2_v1'

     do iNa = 1, nNa
        LogNa = LogNaMin + (LogNaMax-LogNaMin)*real(iNa-1)/real(nNa-1)
        NaSi  = 10.0**LogNa
        Rho   = NaSi * Amean * cAtomicMass

        do iT = 1, nT
           T_K  = TMin * exp(dLogT*real(iT-1))
           ! Use the same exact SI Boltzmann constant as Chromosphere2026.
           ! CRASH's legacy cKToEV is based on a rounded 1.3807e-23 and causes
           ! a resolvable Saha/energy mismatch through exp(-chi/kT).
           TeEV = T_K * BoltzmannSI/cEV

           ! first point of every density column: no previous-Te guess
           UsePreviousTe = (iT > 1)

           call set_ionization_equilibrium(TeEV, NaSi, iError)
           if(iError /= 0)then
              write(*,'(a,i0,a,es14.6,a,es14.6)') &
                   'CRASH EOS error ', iError, ' at T=', T_K, ' n_H=', NaSi
              error stop 1
           end if

           call get_gamma(GammaOut=Gamma,  GammaSOut=GammaS, &
                          GammaeOut=Gammae, GammaSeOut=GammaSe)

           P       = pressure()
           Edens   = NaSi * cEV * internal_energy()
           Cv      = heat_capacity()
           ZbarOut = zAv

           write(iUnit,'(11(1x,es24.16))') &
                T_K, NaSi, Rho, ZbarOut, P, Edens, Gamma, GammaS, Gammae, GammaSe, Cv
           write(iRuntime,'(3(1x,es23.15))') log10(T_K), log10(NaSi), GammaS
        end do
        write(iUnit,'(a)') ''        ! blank line between density blocks (gnuplot friendly)
        write(iRuntime,'(a)') ''
     end do

     close(iUnit)
     close(iRuntime)
     write(*,'(a)') 'wrote '//trim(FileValidation)
     write(*,'(a)') 'wrote '//trim(FileRuntime)

  write(*,'(a)') 'done.'

end program tabulate_gamma
