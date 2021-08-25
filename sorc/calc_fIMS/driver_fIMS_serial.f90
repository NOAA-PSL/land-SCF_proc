 program driver_fIMS

! calculate fractional snow cover on model grid, from IMS snow cover obs. 
! then calculate SWE from fractional snow cover, assuming the snow 
! depletion curve used by the noah model. 
! 
! Clara Draper, July 2021 (based on code from Tseganeh Gichamo, Youlong Xia)


 use IMSaggregate_mod

 implicit none

 integer             :: idim, jdim
 character(len=10)   :: date_str ! yyyymmddhh
 character(len=500)  :: IMS_snowcover_path, IMS_index_path
 logical             :: file_exists
 integer             :: io, ierr

 namelist/fIMS_nml/  idim, jdim, date_str, IMS_snowcover_path, IMS_index_path

 ! read namelist
 inquire(file='fims.nml', exist=file_exists)

 if (.not. file_exists) then
        print *, 'namelistfile does not exist, exiting' 
        stop
 endif

 open (action='read', file='fims.nml', iostat=ierr, newunit=io)
 read (nml=fIMS_nml, iostat=ierr, unit=io) 
 close (io) 

 call calculate_IMS_fsca(idim, jdim, date_str, IMS_snowcover_path, IMS_index_path)

 end program driver_fIMS
