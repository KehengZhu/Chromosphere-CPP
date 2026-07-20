!  Tabulate the effective adiabatic / polytropic index gamma of a partially
!  ionized plasma using the CRASH statistical-sum equation of state
!  (SWMF/util/CRASH).  Written for the Chromosphere2026 project.
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
!  LTE material physics switched on (matching util/CRASH/src/save_eos_table.f90):
!     - ideal ion + electron translation
!     - ionization energy (Saha)          -> the big gamma dip
!     - bound-state excitation (H, He tabulated)   UseExcitation
!     - electron Fermi degeneracy                  UseFermiGas
!     - Coulomb / Debye correction                 UseCoulombCorrection
!
!  Not included in gamma: radiative cooling, radiation energy, conduction,
!  photoionization, or finite-rate/non-equilibrium ionization and recombination.
!  Do not combine this equilibrium-ionization energy closure unchanged with a
!  separate source stage that also pays/returns the ionization potential.
!
!  Two materials are written: pure hydrogen, and a H(0.9)/He(0.1) number mix.

program tabulate_gamma

  use CRASH_ModStatSum
  use CRASH_ModPartition,  ONLY: UseCoulombCorrection, set_mixture, zAv
  use CRASH_ModExcitationData, ONLY: UseExcitation
  use CRASH_ModFermiGas,   ONLY: UseFermiGas, LogGeMinFermi
  use CRASH_ModAtomicMass, ONLY: cAtomicMass_I
  use ModConst,            ONLY: cKToEV, cEV, cAtomicMass

  implicit none

  ! ---- grid ---------------------------------------------------------------
  integer, parameter :: nT  = 360     ! temperature points  (log spaced)
  integer, parameter :: nNa = 46      ! density points      (log spaced)
  real,    parameter :: TMin = 2.0e3, TMax = 5.0e6           ! [K]
  real,    parameter :: LogNaMin = 15.0, LogNaMax = 24.0     ! log10(Na [m^-3])

  ! ---- material description ----------------------------------------------
  integer, parameter :: nMat = 2
  character(len=16)  :: NameMat_I(nMat)  = [ 'hydrogen        ', 'H90He10         ' ]
  character(len=64)  :: FileMat_I(nMat)  = [ &
       'outputs/eos_gamma/gamma_hydrogen.dat                            ', &
       'outputs/eos_gamma/gamma_H90He10.dat                             ' ]

  integer :: iMat, iT, iNa, iError, iUnit
  real    :: T_K, TeEV, NaSi, LogNa, Rho, Amean
  real    :: P, Edens, Gamma, GammaS, Gammae, GammaSe, Cv, ZbarOut
  real    :: dLogT

  !-------------------------------------------------------------------------
  ! full-physics EOS (same switches as save_eos_table.f90)
  UseExcitation        = .true.
  UseFermiGas          = .true.
  UseCoulombCorrection = .true.
  LogGeMinFermi        = -4.0
  UsePreviousTe        = .true.

  dLogT = log(TMax/TMin) / real(nT-1)

  iUnit = 20
  do iMat = 1, nMat

     ! --- set up the material and its mean atomic weight -------------------
     if(iMat == 1)then
        call set_element(1)                                  ! pure hydrogen
        Amean = cAtomicMass_I(1)
     else
        call set_mixture(2, [1,2], [0.9, 0.1])               ! H 90% / He 10%
        Amean = 0.9*cAtomicMass_I(1) + 0.1*cAtomicMass_I(2)
     end if

     open(iUnit, file=trim(FileMat_I(iMat)), status='replace')
     write(iUnit,'(a)') '# CRASH statistical-sum EOS  --  effective gamma table'
     write(iUnit,'(a)') '# material = '//trim(NameMat_I(iMat))// &
          '   (excitation+Fermi+Coulomb ON)'
     write(iUnit,'(a)') '# Gamma  = 1 + P/Edens        (energy gamma, closes e = P/(Gamma-1))'
     write(iUnit,'(a)') '# GammaS = (dlnP/dlnrho)_S     (adiabatic / sound-speed gamma)'
     write(iUnit,'(a)') '# columns:'
     write(iUnit,'(a)') '#   1 T[K]  2 Na[m^-3]  3 Rho[kg/m^3]  4 Zbar  5 P[Pa]  '// &
          '6 Edens[J/m^3]  7 Gamma  8 GammaS  9 Gammae  10 GammaSe  11 Cv[k_B/atom]'

     do iNa = 1, nNa
        LogNa = LogNaMin + (LogNaMax-LogNaMin)*real(iNa-1)/real(nNa-1)
        NaSi  = 10.0**LogNa
        Rho   = NaSi * Amean * cAtomicMass

        do iT = 1, nT
           T_K  = TMin * exp(dLogT*real(iT-1))
           TeEV = T_K * cKToEV

           ! first point of every density column: no previous-Te guess
           UsePreviousTe = (iT > 1)

           call set_ionization_equilibrium(TeEV, NaSi, iError)

           call get_gamma(GammaOut=Gamma,  GammaSOut=GammaS, &
                          GammaeOut=Gammae, GammaSeOut=GammaSe)

           P       = pressure()
           Edens   = NaSi * cEV * internal_energy()
           Cv      = heat_capacity()
           ZbarOut = zAv

           write(iUnit,'(11(1x,es14.6))') &
                T_K, NaSi, Rho, ZbarOut, P, Edens, Gamma, GammaS, Gammae, GammaSe, Cv
        end do
        write(iUnit,'(a)') ''        ! blank line between density blocks (gnuplot friendly)
     end do

     close(iUnit)
     write(*,'(a)') 'wrote '//trim(FileMat_I(iMat))
  end do

  write(*,'(a)') 'done.'

end program tabulate_gamma
