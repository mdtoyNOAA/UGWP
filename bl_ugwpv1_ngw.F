module bl_ugwpv1_ngw

! This module contains the UGWPv1 non-stationary gravity wave drag (NGW) scheme.
!      1) The "V1 CIRES UGWP" scheme as tested in the FV3GFSv16-127L atmosphere model and workflow, which includes
!               the v1 CIRES ugwp non-stationary GW scheme, new revision that generate realistic climate of FV3GFS-127L
!               in the strato-mesosphere in the multi-year simulations (Annual cycles, SAO and QBO in th tropical dynamics).
! See Valery Yudin's presentation at 2020 UFS User's meeting (Jul 2020):
! Gravity waves (GWs): Mesoscale GWs transport momentum, energy (heat) , and create eddy mixing in the whole atmosphere domain; Breaking and dissipating GWs deposit: (a) momentum; (b) heat (energy); and create (c) turbulent mixing of momentum, heat, and tracers
! To properly incorporate GW effects (a-c) unresolved by DYCOREs we need GW physics
! "Unified": a) all GW effects due to both dissipation/breaking; b) identical GW solvers for all GW sources; c) ability to replace solvers.
! Unified Formalism:
! 1. GW Sources: Stochastic and physics based mechanisms for GW-excitations in the lower atmosphere, calibrated by the high-res analyses/forecasts, and observations (3 types of GW sources: orography, convection, fronts/jets).
! 2. GW Propagation: Unified solver for "propagation, dissipation and breaking" excited from all type of GW sources.
! 3. GW Effects: Unified representation of GW impacts on the "resolved" flow for all sources (energy-balanced schemes for momentum, heat and mixing).
! https://www.weather.gov/media/sti/nggps/Presentations%202017/02%20NGGPS_VYUDIN_2017_.pdf
!
!


use mpas_log, only : mpas_log_write
use ccpp_kinds, only: kind_phys
use cires_ugwpv1_module,   only:  cires_ugwpv1_init, ngwflux_update



contains


subroutine ugwpv1_ngw_init (lat_r,levs,dtp,rdzw,dzu,jindx1_tau,jindx2_tau,ddy_j1tau,ddy_j2tau) 

   use ugwp_common
   use mpas_atmphys_constants, only : P0
   use mpas_constants, only : pii, gravity, omega, a, cp, rgas, rv_moist=>rv
   use cires_tauamf_data,   only:  cires_indx_ugwp,  read_tau_amf
   use cires_tauamf_data,   only:  tau_limb, days_limb, ugwp_taulat

   implicit none

   real(kind=kind_phys), dimension(:), intent(in) :: lat_r  ! latitude in radians
   integer, intent(in) :: levs     ! number of model levels
   real(kind=kind_phys), intent(in) :: dtp     ! physics time step
   real(kind=kind_phys), dimension(:), intent(in) :: rdzw ! inverse delta-zeta at u-levels
   real(kind=kind_phys), dimension(:), intent(in) :: dzu  ! delta-zeta at w-levels
   integer, dimension(:), intent(inout) :: jindx1_tau, jindx2_tau
   real(kind=kind_phys), dimension(:), intent(inout) :: ddy_j1tau, ddy_j2tau

   integer :: ios
   logical :: exists
   real    :: dxsg
   integer :: k
   ! Temporary line
   integer :: i

   real(kind=kind_phys) :: p0_std

   real(kind=kind_phys), dimension(levs) :: zu  ! vert. coord. zeta at u-levels centers (m)

   call mpas_log_write('Ahoj svete:  Here we are in ugwpv1_ngw_init!',masterOnly=.true.)

   pi    = pii
   arad  = a        ! Radius of earth (m)
   p0_std = 101325._kind_phys  ! Standard atmospheric pressure (Pa) - Note P0 = ref pressure (10e5)
   grav  = gravity
   omega1= omega
   cpd   = cp
   rd    = rgas
   rv    = rv_moist
   fv    = rv/rd-1._kind_phys

   grav2  = grav + grav; rgrav  = 1._kind_phys/grav ; rgrav2 = rgrav*rgrav
   rdi    = 1._kind_phys / rd ; rcpd = 1._kind_phys/cpd
   rcpd2  = 0.5_kind_phys/cpd
   gor    = grav/rd
   gr2    = grav*gor
   grcp   = grav*rcpd
   gocp   = grcp
   rcpdl  = cpd*rgrav
   grav2cpd = grav*grcp

   pi2      = 2._kind_phys*pi ;  pih = .5_kind_phys*pi
   rad_to_deg=180.0_kind_phys/pi
   deg_to_rad=pi/180.0_kind_phys

   bnv2min = (pi2/1800._kind_phys)*(pi2/1800._kind_phys)
   bnv2max = (pi2/30._kind_phys)*(pi2/30._kind_phys)
   dw2min  = 1.0_kind_phys
   velmin  = sqrt(dw2min)
   minvel  = 0.5_kind_phys

   omega2  = 2._kind_phys*omega1
   omega3  = 3._kind_phys*omega1

   hpscale = 7000._kind_phys ; hpskm = hpscale*1.e-3_kind_phys
   rhp     = 1._kind_phys/hpscale
   rhp2 = 0.5_kind_phys*rhp; rh4 = 0.25_kind_phys*rhp
   rhp4 = rhp2 * rhp2
   khp  = rhp* rd/cpd
   mkzmin  = pi2/80.0e3_kind_phys
   mkz2min = mkzmin*mkzmin
   mkzmax  = pi2/500._kind_phys
   mkz2max = mkzmax*mkzmax
   cdmin   = 2.e-2_kind_phys/mkzmax

   rcpdt  = rcpd/dtp

   ! Create array of vertical coordinate zeta values at layer centers
   zu(1) = 0.5_kind_phys/rdzw(1)     
   do k = 2,levs
      zu(k) = zu(k-1) + dzu(k)
   enddo

   call cires_ugwpv1_init (levs, zu, p0_std, dtp)

   ! Read in tau_amf (momentum flux) for NGW
   call read_tau_amf()

   ! Calculate spatial interpolation weights for NGW tau_amf
   call cires_indx_ugwp(pi,lat_r,jindx1_tau,jindx2_tau,ddy_j1tau,ddy_j2tau)

   ! Temporary lines
   do k = 1,size(lat_r)
      call mpas_log_write('Ahoj svete:   $r $i $i $r $r',      &
              realArgs=(/lat_r(k),ddy_j1tau(k),ddy_j2tau(k)/), &
              intArgs=(/jindx1_tau(k),jindx2_tau(k)/))
   enddo



end subroutine ugwpv1_ngw_init





end module bl_ugwpv1_ngw
