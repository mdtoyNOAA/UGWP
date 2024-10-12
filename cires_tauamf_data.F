module cires_tauamf_data

use ccpp_kinds, only: kind_phys
!...........................................................................................
! tabulated GW-sources: GRACILE/Ern et al., 2018 and/or Resolved GWs from C384-Annual run
!...........................................................................................
implicit none

   integer           :: ntau_d1y, ntau_d2t  
   real(kind=kind_phys), allocatable :: ugwp_taulat(:)
   real(kind=kind_phys), allocatable :: tau_limb(:,:), days_limb(:)
   logical           :: flag_alloctau = .false.          
   character(len=255):: ugwp_taufile =  'ugwp_limb_tau.nc' 

   public :: read_tau_amf, cires_indx_ugwp, tau_amf_interp 

contains
  
    subroutine read_tau_amf()
  
    use  netcdf
    use  mpas_log,           only : mpas_log_write
    use  mpas_derived_types, only : MPAS_LOG_ERR, MPAS_LOG_CRIT
    integer :: ncid,  iernc, vid, dimid, status         
    integer :: k
    
!
      
      iernc=NF90_OPEN(trim(ugwp_taufile), nf90_nowrite, ncid)
     
      if(iernc.ne.0) then         
          call mpas_log_write('read_tau_amf: cannot open file_limb_tab &
                               &data-file '//trim(ugwp_taufile), messageType=MPAS_LOG_CRIT) 
      else
         status = nf90_inq_dimid(ncid, "lat", DimID)
         status = nf90_inquire_dimension(ncid, DimID,  len =ntau_d1y )
         status = nf90_inq_dimid(ncid, "days", DimID)
         status = nf90_inquire_dimension(ncid, DimID,  len =ntau_d2t )
       
         call mpas_log_write('Dimensions of tau_ngw ugwp-v1 $i $i',           &
                              intArgs=(/ntau_d1y,ntau_d2t/),masterOnly=.true.)
         if (ntau_d2t .le. 0 .or. ntau_d1y .le. 0) then 
            call mpas_log_write('ugwp-v1 tau-file='//trim(ugwp_taufile))
            call mpas_log_write('  ugwp-v1:  ntau_d2t= $i, ntau_d1y=$i',      &
                     intArgs=(/ntau_d2t,ntau_d1y/), messageType=MPAS_LOG_CRIT)
         endif
 
         if (.not.allocated(ugwp_taulat))  allocate (ugwp_taulat(ntau_d1y ))
         if (.not.allocated(days_limb))    allocate (days_limb(ntau_d2t))
         if (.not.allocated(tau_limb))     allocate (tau_limb(ntau_d1y, ntau_d2t )) 
              
         iernc = nf90_inq_varid( ncid, 'DAYS', vid )
         iernc = nf90_get_var( ncid, vid, days_limb)
         iernc = nf90_inq_varid( ncid, 'LATS', vid )
         iernc = nf90_get_var( ncid, vid, ugwp_taulat)
         iernc = nf90_inq_varid( ncid, 'ABSMF', vid )
         iernc = nf90_get_var( ncid, vid, tau_limb)
         iernc = nf90_close(ncid)
      endif    

    end subroutine read_tau_amf  

  
    subroutine cires_indx_ugwp (pi,lat_r,j1_tau,j2_tau,w1_j1tau,w2_j2tau)
     
    use ccpp_kinds, only: kind_phys

    implicit none
    
      real(kind=kind_phys) ,   intent(in)                   ::  pi
      real(kind=kind_phys) ,   dimension(:), intent(in)     ::  lat_r    ! latitude in radians
           
      integer, dimension(:), intent(inout)                  ::  j1_tau,   j2_tau
      real(kind=kind_phys) ,   dimension(:), intent(inout)  ::  w1_j1tau, w2_j2tau
      
!locals

      integer :: i,j, j1, j2     
      real(kind=kind_phys), dimension(size(lat_r)) :: dlat    ! latitude in degrees
      real(kind=kind_phys) :: rad2deg
!     
      rad2deg = 180._kind_phys/pi
      dlat(:) = rad2deg*lat_r(:)    ! Calculate latitude in degrees
      do j=1,size(lat_r)
        j2_tau(j) = ntau_d1y
        do i=1,ntau_d1y
          if (dlat(j) < ugwp_taulat(i)) then
            j2_tau(j) = i
            exit
          endif
        enddo
      
        j2_tau(j) = min(j2_tau(j),ntau_d1y)
        j1_tau(j) = max(j2_tau(j)-1,1)

        if (j1_tau(j) /= j2_tau(j) ) then
          w2_j2tau(j) = (dlat(j)  - ugwp_taulat(j1_tau(j))) &
                 / (ugwp_taulat(j2_tau(j))-ugwp_taulat(j1_tau(j)))
        else
          w2_j2tau(j) = 1.0
        endif
          w1_j1tau(j) = 1.0 - w2_j2tau(j)
      enddo
      return
    end subroutine cires_indx_ugwp   

    
    subroutine tau_amf_interp(im,fddd,j1_tau,j2_tau,ddy_j1,ddy_j2,tau_ddd)    
    use ccpp_kinds, only: kind_phys         
    use  mpas_log,           only : mpas_log_write
    use  mpas_derived_types, only : MPAS_LOG_ERR, MPAS_LOG_CRIT
    implicit none
    
!input    
    integer, intent(in)               :: im
    real(kind=kind_phys), intent(in)  :: fddd
      
    real(kind=kind_phys), intent(in),  dimension(im) ::  ddy_j1, ddy_j2
    integer             , intent(in),  dimension(im) ::  j1_tau,j2_tau        
!ouput    
    real(kind=kind_phys), intent(out), dimension(im) ::  tau_ddd
!locals

    integer :: i, j1, j2, it1, it2 , iday
    real(kind=kind_phys)  :: tx1, tx2, w1, w2
    
         it1 = 2
         do iday=1, ntau_d2t
            if (fddd .lt. days_limb(iday) ) then
               it2 = iday
               exit
            endif
         enddo

         it2 = min(it2,ntau_d2t)
         it1 = max(it2-1,1)
         if (it2 > ntau_d2t ) then
            call mpas_log_write(' Error in time-interpolation for tau_amf_interp ', &
                                  messageType=MPAS_LOG_ERR)
            call mpas_log_write(' it1, it2, ntau_d2t   $i $i $i ',                  &
                                  intArgs=(/it1,it2,ntau_d2t/),messageType=MPAS_LOG_ERR)
            call mpas_log_write(' Error in time-interpolation -- see cires_tauamf_data.F90 ', &
                                  messageType=MPAS_LOG_CRIT)
         endif
 
         w2 = (fddd-days_limb(it1))/(days_limb(it2)-days_limb(it1))
         w1 = 1.0-w2     
       
         do i=1, im
            j1 = j1_tau(i)
            j2 = j2_tau(i)
            tx1 = tau_limb(j1, it1)*ddy_j1(i)+tau_limb(j2, it1)*ddy_j2(i)
            tx2 = tau_limb(j1, it2)*ddy_j1(i)+tau_limb(j2, it2)*ddy_j2(i)
            tau_ddd(i) =  tx1*w1 + w2*tx2
         enddo

             
    end subroutine tau_amf_interp  
    
end  module cires_tauamf_data
