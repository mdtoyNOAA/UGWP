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






contains


   subroutine ugwpv1_ngw_init

   implicit none

   print *, "Ahoj svete:  Here we are in ugwpv1_ngw_init!"



   end subroutine ugwpv1_ngw_init





end module bl_ugwpv1_ngw
