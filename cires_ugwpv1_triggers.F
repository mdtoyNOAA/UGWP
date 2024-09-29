module cires_ugwpv1_triggers

use ccpp_kinds, only: kind_phys

contains

!      
!> modified for FV3GFS-127L/C96 QBO-experiments
!! GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch)
   subroutine slat_geos5_2020(im, tau_amp, xlatdeg, tau_gw)

      implicit none
      integer, intent(in) :: im     
      real(kind=kind_phys), intent(in)   :: tau_amp
      real(kind=kind_phys), intent(in)   :: xlatdeg(:)
      real(kind=kind_phys), intent(out)  :: tau_gw(:)
      real(kind=kind_phys)               :: latdeg, flat_gw, tem
      real(kind=kind_phys), parameter    :: fampqbo = 1.25    ! 1.5     
      real(kind=kind_phys), parameter    :: famp60S = 1.0     ! 1.5
      real(kind=kind_phys), parameter    :: famp60N = 1.0     ! 1.0
      real(kind=kind_phys), parameter    :: famp30  = 0.25    ! 0.4
        
      real(kind=kind_phys), parameter    :: swid15  = 12.5  
      real(kind=kind_phys), parameter    :: swid60S = 30.0    ! 40
      real(kind=kind_phys), parameter    :: swid60N = 25.0    ! 30                 
      integer :: i
!    
!
!
      do i=1, im
      
        latdeg = abs(xlatdeg(i))    
        if (latdeg < 15.3) then
          tem = (latdeg-3.0) / swid15
          flat_gw = fampqbo * exp(-tem * tem)
          if (latdeg <= 3.0) flat_gw = fampqbo
        elseif (latdeg <  31.0 .and. latdeg >=  15.3) then
           flat_gw =  famp30
        elseif (latdeg <  60.0 .and. latdeg >=  31.0) then
          tem = (latdeg-60.0) / 23.0
          flat_gw =   famp60N* exp(- tem * tem)
        elseif (latdeg >=  60.0) then
          tem = (latdeg-60.0) /swid60N
          flat_gw =  famp60N * exp(- tem * tem)
        endif

          if (xlatdeg(i) <= -31.0) then   
!                 
            if (latdeg <  60.0 .and. latdeg >=  31.0) then
              tem = (latdeg-60.0) / 23.0
              flat_gw =  famp60S * exp(- tem * tem)           
            endif           
            if (latdeg >=  60.0) then
               tem = (latdeg-60.0) /swid60S
               flat_gw =  famp60S * exp(- tem * tem)
            endif
            
          endif
        tau_gw(i) = tau_amp*flat_gw 
      enddo
!      
   end subroutine slat_geos5_2020   
      

end module cires_ugwpv1_triggers
