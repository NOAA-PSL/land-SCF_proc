 program driver_fSCF

! calculate fractional snow cover on model grid, from IMS or VIIRS snow cover obs. 
! then calculate SWE from fractional snow cover, assuming the snow 
! depletion curve used by the noah model. 
! 
! Clara Draper, July 2021 (based on code from Tseganeh Gichamo, Youlong Xia)
! Yuan Xue, July 2024 (merge VIIRS scf calculation capability with IMS')

 use aggregate_mod, only: calculate_scf

 implicit none

 integer             :: idim, jdim, source
 character(len=11)   :: yyyymmddhh
 character(len=7)    :: jdate
 character(len=20)   :: otype ! orography type, format C$RES (atm) or C${RES}.mx100 (coupled atm/ocean)
 character(len=200)  :: IMS_obs_path, IMS_ind_path, fcst_path
 character(len=200)  :: VIIRS_obs_path, VIIRS_ind_path
 logical             :: file_exists, skip_SD, IS_IMS, IS_VIIRS
 integer             :: io, ierr, lsm, imsformat
 character(len=10)   :: imsversion, imsres, viirsversion
 real                :: viirs_threshold

 namelist/fSCF_nml/  idim, jdim, source, otype, yyyymmddhh, jdate, fcst_path, lsm, skip_SD, & !required input
         IMS_obs_path, IMS_ind_path, imsformat, imsversion, imsres, & !(optional) input needed for IMS DA
         VIIRS_obs_path, VIIRS_ind_path, viirsversion, viirs_threshold !(optional) input needed for VIIRS DA

 ! default values
 source=0
 fcst_path="./"
 lsm=2 ! 1 - noah, 2 - noah-MP. 
 skip_SD=.false.
 IS_IMS=.false.
 IS_VIIRS=.false.

 ! read namelist to get observed snow cover fraction source
 inquire(file='fscf.nml', exist=file_exists)
 if (file_exists) then

        print *, 'Read SCF namelistfile to get the observation source '
        open (action='read', file='fscf.nml', iostat=ierr, newunit=io)
        read (nml=fSCF_nml, iostat=ierr, unit=io)
        close (io) 
        
        select case (source)
        case (1)
            print *, 'Read IMS related inputs ...'
            IS_IMS=.true.
        case (2)
            print *, 'Read VIIRS related inputs ...'
            IS_VIIRS=.true.
        case default
            print *, 'Unrecognized integer for source inputs, exiting ...'
            stop 20
        end select
 else
        print *, 'SCF namelistfile does not exist, exiting ... '
        stop 10
 endif


 if (IS_IMS) then

        call calculate_scf(idim, jdim, source, otype,yyyymmddhh, jdate, fcst_path, lsm, skip_SD, &
                           IMS_obs_path=IMS_obs_path, IMS_ind_path=IMS_ind_path, &
                           imsformat=imsformat, imsversion=imsversion, imsres=imsres)
 endif

 if (IS_VIIRS) then

        call calculate_scf(idim, jdim, source, otype,yyyymmddhh,jdate, fcst_path, lsm, skip_SD, &
                           VIIRS_obs_path=VIIRS_obs_path, VIIRS_ind_path=VIIRS_ind_path, &
                           viirsversion=viirsversion, viirs_threshold=viirs_threshold)
 endif


 end program driver_fSCF
