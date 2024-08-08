module SCFaggregate_mod

use netcdf
use hdf5

private
public calculate_scf

real, parameter    ::  nodata_real = -999. 
integer, parameter ::  nodata_int = -999
real, parameter    ::  nodata_tol = 0.1

! IMS/VIIRS noah-MP snow depth retrieval parameters
real, parameter :: trunc_scf = 0.95 ! For the Noah-MP snow depletion curve, SCF asymptotes to 1. as SD increases 
                                    ! use this value when calculating SD to represent "full" coverage
real, parameter :: snd_max = 300. ! maximum snow depth for snow depth derived from SCF.


! snow depletion curve parameters for IGBP snow depletion curve.
! ORIGINAL
!real, dimension(20), parameter ::  & 
!    mfsno_table = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 2.00, 2.00, &
!                     2.00, 2.00, 2.00, 3.00, 3.00, 4.00, 4.00, &
!                     2.50, 3.00, 3.00, 3.50, 3.50, 3.50 /)

! LIMIT MFSNO to 3.
real, dimension(20), parameter ::  & 
    mfsno_table = (/ 1.00, 1.00, 1.00, 1.00, 1.00, 2.00, 2.00, &
                     2.00, 2.00, 2.00, 3.00, 3.00, 3.00, 3.00, &
                     2.50, 3.00, 3.00, 3.00, 3.00, 3.00 /)

real, dimension(20), parameter ::  & 
    scffac_table = (/ 0.005, 0.005, 0.005, 0.005, 0.005, 0.008, &
                      0.008, 0.010, 0.010, 0.010, 0.010, 0.007, 0.021, &  
                      0.013, 0.015, 0.008, 0.015, 0.015, 0.015, 0.015 /)

contains

!====================================
! main routine to read in inputs, calculate IMS/VIIRS snow cover fraction, IMS/VIIRS SWE, 
! then IMS/VIIRS SND, and write out results on model grid.
! SWE is calculated using the model relationship. 
! SD is calculated using the forecast snow density. 
! SD is QC'ed out where both model and obs have "full" snow cover 
! (since can get no info from IMS/VIIRS snow cover in this case)

! note: original version of this code used model and observed snow cover as a fraction. 
! However, comparison between model and observations shows that the obs can have <100% 
! snow cover for very deep snow (likely, due to satellite seeing trees, etc), while this 
! does not (cannot) occur in the model.  This resulted in the DA incorrectly removing snow 
! when this occured. The code has now been reverted to process the IMS/VIIRS SCF as 100% for SFC>50%, and 
! as 0% for SCF < 50%.

! note on timing: files are once a day, nominally at 0 UTC. 
subroutine calculate_scf(idim, jdim, source, otype, yyyymmddhh, jdate, fcst_path, lsm, skip_SD, &
                         IMS_obs_path, VIIRS_obs_path, &
                         IMS_ind_path, VIIRS_ind_path, imsformat, imsversion, viirsversion, viirs_threshold, imsres) 
        implicit none
        
        integer, intent(in)                      :: idim, jdim, lsm, source !1=IMS; 2=VIIRS
        integer, intent(in),optional             :: imsformat
        character(len=20), intent(in)            :: otype  
        character(len=11), intent(in)            :: yyyymmddhh
        character(len=7), intent(in)             :: jdate
        character(len=*), intent(in),optional    :: IMS_obs_path, IMS_ind_path
        character(len=*), intent(in),optional    :: VIIRS_obs_path, VIIRS_ind_path
        character(len=*), intent(in)             :: fcst_path
        logical, intent(in)                      :: skip_SD
        character(len=10), intent(in),optional   :: imsversion, viirsversion, imsres 
        real, optional                           :: viirs_threshold

        real                :: vtype(idim,jdim,6)       ! model vegetation type
        integer             :: landmask(idim,jdim,6)
        real                :: swefcs(idim,jdim,6), sndfcs(idim,jdim,6) ! forecast SWE, SND
        real                :: stcfcs(idim,jdim,6) ! forecast soil temp  
        real                :: denfcs(idim,jdim,6) ! forecast density
        real                :: scffcs(idim,jdim,6) ! forecast snow cover fraction
        real                :: scf(idim,jdim,6) ! IMS/VIIRS snow cover fraction, on model grid
        real                :: swe(idim,jdim,6) ! SWE derived from scf, on model grid - only needed for noah lsm
        real                :: snd(idim,jdim,6) ! snow depth derived from scf, on model grid
        real                :: lonFV3(idim,jdim,6) ! longitude on model grid
        real                :: latFV3(idim,jdim,6) ! latutide on model grid
        real                :: oroFV3(idim,jdim,6) ! orography model grid
        character(len=250)  :: IMS_obs_file
        character(len=250)  :: VIIRS_obs_file
        character(len=8)    :: date_from_file
        integer             :: i,j,t
        integer             :: time
!=============================================================================================
! 1. Read forecast info, and IMS/VIIRS data and indexes from file, then calculate SWE
!=============================================================================================

        if (.not. skip_SD) then
            ! note: snow cover being read here is calculated at the start of the 
            !       time step, and does not account for snow changes over the last 
            !       time step. Resulting errors are generally very small.
            call  read_fcst(fcst_path, yyyymmddhh, idim, jdim, vtype, swefcs,  & 
                            sndfcs, stcfcs, landmask)

            call calc_density(idim, jdim, lsm, landmask, swefcs, sndfcs, stcfcs, denfcs)
        endif 

        ! read in either ascii or nc IMS obs, and indexes, map to model grid

        if (source .eq. 1) then
          if (imsformat==1) then
            IMS_obs_file = trim(IMS_obs_path)//"ims"//trim(jdate)//"_"//trim(imsres)//"_v"//trim(imsversion)//".asc"
          elseif (imsformat==2) then  
            IMS_obs_file = trim(IMS_obs_path)//"ims"//trim(jdate)//"_"//trim(imsres)//"_v"//trim(imsversion)//".nc"
          else
            print *, 'fatal error reading IMS snow cover data '    
          endif
       
          print *, 'reading IMS snow cover data from ', trim(IMS_obs_file) 

          call read_IMS_onto_model_grid(IMS_obs_file, IMS_ind_path, imsformat, imsres, &
                                   jdim, idim, otype, lonFV3, latFV3, oroFV3, scf, date_from_file, time)
        endif

        ! read in hdf5 VIIRS obs, and indexes, map to model grid

        if (source .eq. 2) then

          VIIRS_obs_file = trim(VIIRS_obs_path)//"VNP10C1.A"//trim(jdate)//"."//trim(viirsversion)//".h5"

          print *, 'reading VIIRS snow cover data from ', trim(VIIRS_obs_file)

          call read_VIIRS_onto_model_grid(VIIRS_obs_file, VIIRS_ind_path, viirs_threshold, &
                                   jdim, idim, otype, vtype, lonFV3, latFV3, oroFV3, scf)
        endif


       if (.not. skip_SD) then
            ! calculate SWE from IMS/VIIRS snow cover fraction (using model relationship)
            ! no value is calculated if both IMS/VIIRS and model have 100% snow cover
            ! also removes scf if model is non-land
            if (lsm==1) then
                call calcSWE_noah(scf, vtype, swefcs, idim, jdim, swe)
                ! calculate snow depth from IMS/VIIRS SWE, using model density
                do t=1,6
                  do i=1,idim
                    do j=1,jdim 
                      if  ( abs( swe(i,j,t) -nodata_real ) > nodata_tol ) then
                        snd(i,j,t) = swe(i,j,t)/denfcs(i,j,t)
                      else
                        snd(i,j,t) = nodata_real
                      endif
                    enddo
                  enddo
                enddo
            elseif (lsm==2) then 
                ! calculate SD from IMS/VIIRS SCF
                call calcSD_noahmp(scf, vtype, denfcs, idim, jdim, snd)

                ! calculate SCF from model forecast SD and SWE (since not always in restart)
                call calcSCF_noahmp(vtype, denfcs, sndfcs, idim, jdim, scffcs) 

                !exclude IMS/VIIRS snow depth, where both IMS/VIIRS and model are 100% 
                do t=1,6
                  do i=1,idim
                    do j=1,jdim 
                            ! do not assimilate, if obs and model indicate "full" snow. 
                            ! (recall: converting obs SCF>0.5 to full snow)
                            ! additional check to remove obs if obs have full snow, but calculated snow depth is lower than the model
                            if ( (scf(i,j,t) >= 0.5 ) .and.  & 
                                    (  (scffcs(i,j,t) > trunc_scf ) .or.  (sndfcs(i,j,t ) >  snd(i,j,t) ) ) ) then 
                                   snd(i,j,t) = nodata_real
                            endif 
                            if ( snd(i,j,t) > snd_max ) then 
                                   snd(i,j,t) = nodata_real
                            endif
                    enddo 
                  enddo 
                enddo

            else 
                print *, 'unknown lsm:', lsm, ', choose 1 - noah, 2 - noah-mp' 
                stop 10 
            endif
        else 
                snd=nodata_real
        endif ! skip_SD
!=============================================================================================
! 2.  Write outputs
!=============================================================================================
       
        if (source .eq. 1) then
           !call write_IMS_outputs_2D(idim, jdim, scf,snd)
           call write_IMS_outputs_vec(idim, jdim, otype, yyyymmddhh, scf, snd, lonFV3, latFV3, oroFV3, date_from_file, time)
        endif

        if (source .eq. 2) then
               call write_VIIRS_outputs_vec(idim, jdim, otype, yyyymmddhh, scf, snd, lonFV3, latFV3, oroFV3)
        endif

        return

 end subroutine calculate_scf

!====================================
! routine to write the output to file - 2D (one file per tile)
! YX: Only available (but unused as of 07/15/2024) for IMS observations

 subroutine write_IMS_outputs_2D(idim, jdim, scf, snd)

      !------------------------------------------------------------------
      !------------------------------------------------------------------
      implicit none

      integer, intent(in)         :: idim, jdim
      real, intent(in)            :: scf(idim,jdim,6)
      real, intent(in)            :: snd(idim,jdim,6)

      character(len=250)          :: output_file
      character(len=1)            :: tile_str
      integer                     :: fsize=65536, inital=0
      integer                     :: header_buffer_val = 16384
      integer                     :: dIMS_3d(3), dIMS_strt(3), dIMS_end(3)
      integer                     :: error, i, ncid
      integer                     :: dim_x, dim_y, dim_time
      integer                     :: id_x, id_y, id_time
      integer                     :: id_scfIMS, id_sndIMS 
      integer                     :: itile
 
      real(kind=4)                :: times
      real(kind=4)                :: xy_data(idim)

      do itile = 1, 6

        write(tile_str, '(i1.1)') itile ! assuming <10 tiles.
        output_file = "./IMSfSCA.tile"//tile_str//".nc"
        print*,'writing output to ',trim(output_file) 
        
        !--- create the file
        error = nf90_create(output_file, ior(nf90_netcdf4,nf90_classic_model), ncid, initialsize=inital, chunksize=fsize)
        call netcdf_err(error, 'creating file='//trim(output_file) )

        !--- define dimensions
        error = nf90_def_dim(ncid, 'xaxis_1', idim, dim_x)
        call netcdf_err(error, 'defining xaxis dimension' )
        error = nf90_def_dim(ncid, 'yaxis_1', jdim, dim_y)
        call netcdf_err(error, 'defining yaxis dimension' )
        error = nf90_def_dim(ncid, 'Time', 1, dim_time)
        call netcdf_err(error, 'defining time dimension' )

        !--- define fields
        error = nf90_def_var(ncid, 'xaxis_1', nf90_float, dim_x, id_x)
        call netcdf_err(error, 'defining xaxis_1 field' )
        error = nf90_put_att(ncid, id_x, "long_name", "xaxis_1")
        call netcdf_err(error, 'defining xaxis_1 long name' )
        error = nf90_put_att(ncid, id_x, "units", "none")
        call netcdf_err(error, 'defining xaxis_1 units' )
        error = nf90_put_att(ncid, id_x, "cartesian_axis", "X")
        call netcdf_err(error, 'writing xaxis_1 field' )

        error = nf90_def_var(ncid, 'yaxis_1', nf90_float, dim_y, id_y)
        call netcdf_err(error, 'defining yaxis_1 field' )
        error = nf90_put_att(ncid, id_y, "long_name", "yaxis_1")
        call netcdf_err(error, 'defining yaxis_1 long name' )
        error = nf90_put_att(ncid, id_y, "units", "none")
        call netcdf_err(error, 'defining yaxis_1 units' )
        error = nf90_put_att(ncid, id_y, "cartesian_axis", "Y")
        call netcdf_err(error, 'writing yaxis_1 field' )

        error = nf90_def_var(ncid, 'Time', nf90_float, dim_time, id_time)
        call netcdf_err(error, 'defining time field' )
        error = nf90_put_att(ncid, id_time, "long_name", "Time")
        call netcdf_err(error, 'defining time long name' )
        error = nf90_put_att(ncid, id_time, "units", "time level")
        call netcdf_err(error, 'defining time units' )
        error = nf90_put_att(ncid, id_time, "cartesian_axis", "T")
        call netcdf_err(error, 'writing time field' )

        dIMS_3d(1) = dim_x
        dIMS_3d(2) = dim_y
        dIMS_3d(3) = dim_time

        error = nf90_def_var(ncid, 'IMSscf', nf90_double, dIMS_3d, id_scfIMS)
        call netcdf_err(error, 'defining IMSscf' )
        error = nf90_put_att(ncid, id_scfIMS, "long_name", "IMS snow covered fraction")
        call netcdf_err(error, 'defining IMSscf long name' )
        error = nf90_put_att(ncid, id_scfIMS, "units", "-")
        call netcdf_err(error, 'defining IMSscf units' )

        error = nf90_def_var(ncid, 'IMSsnd', nf90_double, dIMS_3d, id_sndIMS)
        call netcdf_err(error, 'defining IMSsnd' )
        error = nf90_put_att(ncid, id_sndIMS, "long_name", "IMS snow depth")
        call netcdf_err(error, 'defining IMSsnd long name' )
        error = nf90_put_att(ncid, id_sndIMS, "units", "mm")
        call netcdf_err(error, 'defining IMSsnd units' )

        error = nf90_enddef(ncid, header_buffer_val,4,0,4)
        call netcdf_err(error, 'defining header' )

        do i = 1, idim
        xy_data(i) = float(i)
        enddo
        times = 1.0

        error = nf90_put_var( ncid, id_x, xy_data)
        call netcdf_err(error, 'writing xaxis record' )
        error = nf90_put_var( ncid, id_y, xy_data)
        call netcdf_err(error, 'writing yaxis record' )
        error = nf90_put_var( ncid, id_time, times)
        call netcdf_err(error, 'writing time record' )

        dIMS_strt(1:3) = 1
        dIMS_end(1) = idim
        dIMS_end(2) = jdim
        dIMS_end(3) = 1
        
        error = nf90_put_var(ncid, id_scfIMS, scf(:,:,itile), dIMS_strt, dIMS_end)
        call netcdf_err(error, 'writing IMSscf record')

        error = nf90_put_var(ncid, id_sndIMS, snd(:,:,itile), dIMS_strt, dIMS_end)
        call netcdf_err(error, 'writing IMSsnd record')

        error = nf90_close(ncid)

      end do
    
 end subroutine write_IMS_outputs_2D


!====================================
! routine to write the output to file - vector
! also writes out the model lat/lon for the grid cell that the data have been 
! processed onto.

 subroutine write_IMS_outputs_vec(idim, jdim, otype, date_str,scf, snd, lonFV3, latFV3, oroFV3, date_from_file, time)

    implicit none
    character(len=11), intent(in)  :: date_str
    character(len=20), intent(in)  :: otype
    character(len=8), intent(in)   :: date_from_file
    integer, intent(in)            :: idim, jdim
    integer, intent(in)            :: time
    real, intent(in)            :: scf(idim,jdim,6)
    real, intent(in)            :: snd(idim,jdim,6)
    real, intent(in)            :: latFV3(idim,jdim,6)
    real, intent(in)            :: lonFV3(idim,jdim,6)
    real, intent(in)            :: oroFV3(idim,jdim,6)

    character(len=250)          :: output_file
    character(len=10)           :: time_char
    integer                     :: header_buffer_val = 16384
    integer                     :: dim_time, id_time
    integer                     :: i,j,t,n, nobs
    integer                     :: error, ncid
    integer                     :: id_scfIMS, id_sndIMS , id_obs, id_lon, id_lat, id_oro
    real, allocatable           :: data_vec(:,:)
    real, allocatable           :: coor_vec(:,:)
 
    output_file = "./IMSscf."//date_str(1:8)//"."//trim(adjustl(otype))//".nc"
    print*,'writing output to ',trim(output_file) 
    
    !--- create the file
    !error = nf90_create(output_file, ior(nf90_netcdf4,nf90_classic_model), ncid, initialsize=inital, chunksize=fsize)
    error = nf90_create(output_file, ior(nf90_netcdf4,nf90_classic_model), ncid)
    call netcdf_err(error, 'creating file='//trim(output_file) )

    ! collect obs
    ! note: writing out all scf ovs. snd will be missing for many locations 
    !       since it is removed if both model and ovs have scf==1
    nobs = count (abs(scf -nodata_real) > nodata_tol)
    print *, 'writing out', nobs, ' observations'

    if(date_from_file == "ascifile") then
      time_char =  date_str(1:8)//"00"
    endif

    allocate(data_vec(2,nobs)) 
    allocate(coor_vec(3,nobs)) 

    !--- define spatial dimension
    error = nf90_def_dim(ncid, 'numobs', nobs, id_obs)
    call netcdf_err(error, 'defining obs dimension' )

    ! --- define global attributes and units for valid_time_str and valid_epoch_time
    if(date_from_file == "ascifile") then

        error = nf90_put_att(ncid, NF90_GLOBAL, "valid_time_str", time_char)
        call netcdf_err(error, 'put valid_date_str as global attribute')

        error = nf90_put_att(ncid, NF90_GLOBAL, "valid_time_str_comment", "This is the 00Z reference time. &
          Note that products are nowcasted to be valid specifically at the time given here.")
        call netcdf_err(error, 'defining valid_time_str comment' )

    else   
        error = nf90_put_att(ncid, NF90_GLOBAL, "valid_epoch_time", time)
        call netcdf_err(error, 'put valid_epoch_time as global attribute')

        error = nf90_put_att(ncid, NF90_GLOBAL, "valid_epoch_time_comment", "This is the 00Z reference time. &
          Note that products are nowcasted to be valid specifically at the time given here.")
        call netcdf_err(error, 'defining "valid_epoch_time comment' )

        error = nf90_put_att(ncid, NF90_GLOBAL, "valid_epoch_time_units", "seconds since 1970-01-01T00:00:00Z")
        call netcdf_err(error, 'defining time units' )    
    endif
    
    !--- define longitude 
    error = nf90_def_var(ncid, 'lon', nf90_double, id_obs, id_lon)
    call netcdf_err(error, 'defining lon' )
    error = nf90_put_att(ncid, id_lon, "long_name", "longitude")
    call netcdf_err(error, 'defining lon long name' )

    !--- define latitude
    error = nf90_def_var(ncid, 'lat', nf90_double, id_obs, id_lat)
    call netcdf_err(error, 'defining lat' )
    error = nf90_put_att(ncid, id_lat, "long_name", "latitude")
    call netcdf_err(error, 'defining lat long name' )

    !--- define orography
    error = nf90_def_var(ncid, 'oro', nf90_double, id_obs, id_oro)
    call netcdf_err(error, 'defining oro' )
    error = nf90_put_att(ncid, id_oro, "long_name", "orography")
    call netcdf_err(error, 'defining oro long name' )

    !--- define snow cover
    error = nf90_def_var(ncid, 'IMSscf', nf90_double, id_obs, id_scfIMS)
    call netcdf_err(error, 'defining IMSscf' )
    error = nf90_put_att(ncid, id_scfIMS, "long_name", "IMS snow covered fraction")
    call netcdf_err(error, 'defining IMSscf long name' )
    error = nf90_put_att(ncid, id_scfIMS, "units", "-")
    call netcdf_err(error, 'defining IMSscf units' )

    !--- define snow depth
    error = nf90_def_var(ncid, 'IMSsnd', nf90_double, id_obs, id_sndIMS)
    call netcdf_err(error, 'defining IMSsnd' )
    error = nf90_put_att(ncid, id_sndIMS, "long_name", "IMS snow depth")
    call netcdf_err(error, 'defining IMSsnd long name' )
    error = nf90_put_att(ncid, id_sndIMS, "units", "mm")
    call netcdf_err(error, 'defining IMSsnd units' )

    error = nf90_enddef(ncid)
    call netcdf_err(error, 'defining header' )

    data_vec=nodata_real
    coor_vec=nodata_real

    n=0
    do t=1,6
     do i=1,idim 
      do j=1,jdim
        if (abs(scf(i,j,t) -nodata_real) > nodata_tol) then 
                n=n+1
                data_vec(1,n) = scf(i,j,t)
                data_vec(2,n) = snd(i,j,t)
                coor_vec(1,n) = lonFV3(i,j,t)
                coor_vec(2,n) = latFV3(i,j,t)
                coor_vec(3,n) = oroFV3(i,j,t)
        endif
      enddo 
     enddo 
    enddo
   
    ! --- put lat, lon, data

    error = nf90_put_var(ncid, id_scfIMS, data_vec(1,:))
    call netcdf_err(error, 'writing IMSscf record')

    error = nf90_put_var(ncid, id_sndIMS, data_vec(2,:))
    call netcdf_err(error, 'writing IMSsnd record')

    error = nf90_put_var(ncid, id_lon, coor_vec(1,:))
    call netcdf_err(error, 'writing lon record')

    error = nf90_put_var(ncid, id_lat, coor_vec(2,:))
    call netcdf_err(error, 'writing lat record')

    error = nf90_put_var(ncid, id_oro, coor_vec(3,:))
    call netcdf_err(error, 'writing oro record')
  
    error = nf90_close(ncid)
    deallocate(data_vec)
    deallocate(coor_vec)
    
 end subroutine write_IMS_outputs_vec

!====================================
! routine to write the VIIRS output to file - vector
! also writes out the model lat/lon for the grid cell that the data have been
! processed onto.

 subroutine write_VIIRS_outputs_vec(idim, jdim, otype, date_str,scf, snd, lonFV3, latFV3, oroFV3)

    implicit none
    character(len=11), intent(in)  :: date_str
    character(len=20), intent(in)  :: otype
    integer, intent(in)         :: idim, jdim
    real, intent(in)            :: scf(idim,jdim,6)
    real, intent(in)            :: snd(idim,jdim,6)
    real, intent(in)            :: latFV3(idim,jdim,6)
    real, intent(in)            :: lonFV3(idim,jdim,6)
    real, intent(in)            :: oroFV3(idim,jdim,6)

    character(len=250)          :: output_file
    integer                     :: header_buffer_val = 16384
    integer                     :: i,j,t,n, nobs
    integer                     :: error, ncid
    integer                     :: id_scfVIIRS, id_sndVIIRS , id_obs, id_lon, id_lat, id_oro
    real, allocatable           :: data_vec(:,:)
    real, allocatable           :: coor_vec(:,:)

    output_file = "./VIIRSscf."//date_str(1:8)//"."//trim(adjustl(otype))//".nc"
    print*,'writing output to ',trim(output_file)

    !--- create the file
    !error = nf90_create(output_file, ior(nf90_netcdf4,nf90_classic_model), ncid, initialsize=inital, chunksize=fsize)
    error = nf90_create(output_file, ior(nf90_netcdf4,nf90_classic_model), ncid)
    call netcdf_err(error, 'creating file='//trim(output_file) )

    ! collect obs
    ! note: writing out all scf ovs. snd will be missing for many locations
    !       since it is removed if both model and ovs have scf==1
    nobs = count (abs(scf -nodata_real) > nodata_tol)
    print *, 'writing out', nobs, ' observations'

    allocate(data_vec(2,nobs))
    allocate(coor_vec(3,nobs))

    !--- define dimension
    error = nf90_def_dim(ncid, 'numobs', nobs, id_obs)
    call netcdf_err(error, 'defining obs dimension' )

    !--- define longitude
    error = nf90_def_var(ncid, 'lon', nf90_double, id_obs, id_lon)
    call netcdf_err(error, 'defining lon' )
    error = nf90_put_att(ncid, id_lon, "long_name", "longitude")
    call netcdf_err(error, 'defining lon long name' )

    !--- define latitude
    error = nf90_def_var(ncid, 'lat', nf90_double, id_obs, id_lat)
    call netcdf_err(error, 'defining lat' )
    error = nf90_put_att(ncid, id_lat, "long_name", "latitude")
    call netcdf_err(error, 'defining lat long name' )

    !--- define orography
    error = nf90_def_var(ncid, 'oro', nf90_double, id_obs, id_oro)
    call netcdf_err(error, 'defining oro' )
    error = nf90_put_att(ncid, id_oro, "long_name", "orography")
    call netcdf_err(error, 'defining oro long name' )

    !--- define snow cover
    error = nf90_def_var(ncid, 'VIIRSscf', nf90_double, id_obs, id_scfVIIRS)
    call netcdf_err(error, 'defining VIIRSscf' )
    error = nf90_put_att(ncid, id_scfVIIRS, "long_name", "VIIRS snow covered fraction")
    call netcdf_err(error, 'defining VIIRSscf long name' )
    error = nf90_put_att(ncid, id_scfVIIRS, "units", "-")
    call netcdf_err(error, 'defining VIIRSscf units' )

    !--- define snow depth
    error = nf90_def_var(ncid, 'VIIRSsnd', nf90_double, id_obs, id_sndVIIRS)
    call netcdf_err(error, 'defining VIIRSsnd' )
    error = nf90_put_att(ncid, id_sndVIIRS, "long_name", "VIIRS snow depth")
    call netcdf_err(error, 'defining VIIRSsnd long name' )
    error = nf90_put_att(ncid, id_sndVIIRS, "units", "mm")
    call netcdf_err(error, 'defining VIIRSsnd units' )

    error = nf90_enddef(ncid)
    call netcdf_err(error, 'defining header' )

    data_vec=nodata_real
    coor_vec=nodata_real

    n=0
    do t=1,6
     do i=1,idim
      do j=1,jdim
        if (abs(scf(i,j,t) -nodata_real) > nodata_tol) then
                n=n+1
                data_vec(1,n) = scf(i,j,t)
                data_vec(2,n) = snd(i,j,t)
                coor_vec(1,n) = lonFV3(i,j,t)
                coor_vec(2,n) = latFV3(i,j,t)
                coor_vec(3,n) = oroFV3(i,j,t)
        endif
      enddo
     enddo
    enddo

    error = nf90_put_var(ncid, id_scfVIIRS, data_vec(1,:))
    call netcdf_err(error, 'writing VIIRSscf record')

    error = nf90_put_var(ncid, id_sndVIIRS, data_vec(2,:))
    call netcdf_err(error, 'writing VIIRSsnd record')

    error = nf90_put_var(ncid, id_lon, coor_vec(1,:))
    call netcdf_err(error, 'writing lon record')

    error = nf90_put_var(ncid, id_lat, coor_vec(2,:))
    call netcdf_err(error, 'writing lat record')

    error = nf90_put_var(ncid, id_oro, coor_vec(3,:))
    call netcdf_err(error, 'writing oro record')

    error = nf90_close(ncid)
    deallocate(data_vec)
    deallocate(coor_vec)

 end subroutine write_VIIRS_outputs_vec

!====================================
! read in required forecast fields from a UFS surface restart 

 subroutine read_fcst(path, date_str, idim, jdim, vetfcs, swefcs, sndfcs, stcfcs, landmask)

        implicit none
        character(len=*), intent(in)      :: path
        character(11), intent(in)         :: date_str
        integer, intent(in)               :: idim, jdim
        real, intent(out)                 :: vetfcs(idim,jdim,6), swefcs(idim,jdim,6)
        real, intent(out)                 :: stcfcs(idim,jdim,6)
        real, intent(out)                 :: sndfcs(idim,jdim,6)
        integer, intent(out)              :: landmask(idim,jdim,6)

        integer                   :: error, ncid, i,j, t 
        integer                   :: id_dim, id_var, idim_file
        character                 :: tt
        character(len=300)        :: fcst_file

        real(kind=8)              :: dummy(idim,jdim)
        real(kind=8)              :: dummy3(idim,jdim,4) ! 4 = number of soil layers
        logical                   :: file_exists

        integer, parameter        :: veg_type_landice = 15

        do t =1,6
            write(tt, "(i1)") t
            fcst_file = trim(path)//trim(date_str)// & 
                                "0000.sfc_data.tile"//tt//".nc"

            print *, 'reading model backgroundfile:', trim(fcst_file)

            inquire(file=trim(fcst_file), exist=file_exists)

            if (.not. file_exists) then
                    print *, 'read_fcst error,file does not exist', &
                            trim(fcst_file) , ' exiting'
                    stop 10
            endif

            error=nf90_open(trim(fcst_file), nf90_nowrite,ncid)
            call netcdf_err(error, 'opening file: '//trim(fcst_file) )

            ! check dimension 
            error=nf90_inq_dimid(ncid, 'xaxis_1', id_dim)
            call netcdf_err(error, 'reading xaxis_1' )
            error=nf90_inquire_dimension(ncid,id_dim,len=idim_file)
            call netcdf_err(error, 'reading xaxis_1' )

            if ((idim_file) /= idim) then
                print*,'fatal error reading fcst file: dimensions wrong.'
                stop 10
            endif

            ! vegetation type
            error=nf90_inq_varid(ncid, "vtype", id_var)
            call netcdf_err(error, 'reading vtype id' )
            error=nf90_get_var(ncid, id_var, dummy)
            call netcdf_err(error, 'reading vtype' )
            vetfcs(:,:,t) = dummy

            ! Snow water equivalent
            error=nf90_inq_varid(ncid, "sheleg", id_var)
            call netcdf_err(error, 'reading sheleg id' )
            error=nf90_get_var(ncid, id_var, dummy)
            call netcdf_err(error, 'reading sheleg' )
            swefcs(:,:,t) = dummy

            ! snow depth
            error=nf90_inq_varid(ncid, "snwdph", id_var)
            call netcdf_err(error, 'reading snwdph id' )
            error=nf90_get_var(ncid, id_var, dummy)
            call netcdf_err(error, 'reading snwdph' )
            sndfcs(:,:,t) = dummy

            ! layer 1 soil temperature
            error=nf90_inq_varid(ncid, "stc", id_var)
            call netcdf_err(error, 'reading stc id' )
            error=nf90_get_var(ncid, id_var, dummy3)
            call netcdf_err(error, 'reading stc' )
            stcfcs(:,:,t) = dummy3(:,:,1) 

            ! land mask
            error=nf90_inq_varid(ncid, "slmsk", id_var)
            call netcdf_err(error, 'reading slmsk id' )
            error=nf90_get_var(ncid, id_var, dummy)
            call netcdf_err(error, 'reading slmsk' )

            ! slmsk in file is: 0 - ocean, 1 - land, 2 -seaice
            ! convert to: integer with  0 - glacier or non-land, 1 - non-glacier covered land

            do i = 1, idim 
              do j = 1, jdim
               ! if land, but not land ice, set mask to 1.
               if ( (nint(dummy(i,j)) == 1 ) .and.   &
                    ( nint(vetfcs(i,j,t)) /=  veg_type_landice  )) then
                    landmask(i,j,t) = 1
               else
                    landmask(i,j,t) = 0
               endif
              enddo
            enddo

        error = nf90_close(ncid)

    enddo

 end subroutine read_fcst

!====================================
! read in the IMS observations and associated index file, then 
! aggregate onto the model grid.

 subroutine read_IMS_onto_model_grid(IMS_obs_file, IMS_ind_path, &
            imsformat, imsres, jdim, idim, otype, lonFV3, latFV3, oroFV3,scf, date_from_file, time)
                    
        implicit none
    
        character(len=*), intent(in)   :: IMS_obs_file, IMS_ind_path
        integer, intent(in)            :: jdim, idim, imsformat
        character(len=20), intent(in)  :: otype
        character(len=10), intent(in)  :: imsres
        character(len=8), intent(out)  :: date_from_file
        real, intent(out)              :: scf(jdim,idim,6)     
        real, intent(out)              :: lonFV3(jdim,idim,6)     
        real, intent(out)              :: latFV3(jdim,idim,6)     
        real, intent(out)              :: oroFV3(jdim,idim,6)     
        integer, intent(out)           :: time  
 
        integer, allocatable    :: IMS_flag(:,:)   
        integer, allocatable    :: IMS_index(:,:,:)
        real                    :: land_points(jdim,idim,6), snow_points(jdim,idim,6)
        
        integer                :: error, ncid, id_dim, id_var, n_ind, id_time
        integer                :: i_ims, j_ims, itile, tile, tile_i, tile_j
        logical                :: file_exists
        character(len=250)     :: IMS_ind_file
        character(len=20)      :: datestring

        integer                :: icol, irow

        ! read IMS observations in
        inquire(file=trim(IMS_obs_file), exist=file_exists)

        if (.not. file_exists) then
           print *, 'observation_read_IMS_full error,file does not exist', &
                        trim(IMS_obs_file) , ' exiting'
           stop 10

        endif

        ! to do - better to read these from file? 
        if (trim(imsres) == "4km" ) then 
            i_ims = 6144
            j_ims = 6144
        elseif (trim(imsres) == "24km" ) then 
            i_ims = 1024
            j_ims = 1024
        else 
           print *, 'unrecognised imsres', trim(imsres), ' exiting'
           stop 10
        endif

        allocate(IMS_flag(j_ims, i_ims))   
        
        if (imsformat==1) then
        ! read in ascii IMS data  
            open(10, file=IMS_obs_file, form="formatted", status="old")
          
            do irow = 1, 30
             read(10,*)     ! read ims ascii header
            end do

            do irow = 1, j_ims
               read(10,'(6144i1)') (IMS_flag(icol, irow), icol=1, i_ims)
            end do
 
            date_from_file='ascifile'
            time=-999999           
 
        elseif (imsformat==2) then
        ! read in netCDF IMS data
           error=nf90_open(trim(IMS_obs_file),nf90_nowrite, ncid)
           call netcdf_err(error, 'opening file: '//trim(IMS_obs_file) )

           error=nf90_inq_varid(ncid, 'IMS_Surface_Values', id_var)
           call netcdf_err(error, 'error reading IMS id' )
           
           error=nf90_get_var(ncid, id_var, IMS_flag)
           call netcdf_err(error, 'error reading IMS nc data' )

           error=nf90_inq_varid(ncid, 'time', id_time)
           call netcdf_err(error, 'error reading time id' )

           error=nf90_get_var(ncid, id_time, time)
           call netcdf_err(error, 'error reading time nc data' )
 
           date_from_file='nc_input'
           
           error = nf90_close(ncid)

        else
           print*,'fatal error reading IMS OBS file'
           stop 10   
        endif 

        ! IMS codes: 0 - outside range,
        !          : 1 - sea
        !          : 2 - land, no snow 
        !          : 3 - sea ice 
        !          : 4 - snow covered land


        where(IMS_flag == 0 ) IMS_flag = nodata_int ! set outside range to NA
        where(IMS_flag == 1 ) IMS_flag = nodata_int ! set sea to NA
        where(IMS_flag == 3 ) IMS_flag = nodata_int ! set sea ice to NA
        where(IMS_flag == 2 ) IMS_flag = 0          ! set land, no snow to 0
        where(IMS_flag == 4 ) IMS_flag = 1          ! set snow on land to 1

        ! read index file for mapping IMS to model grid 

        IMS_ind_file = trim(IMS_ind_path)//"IMS"//trim(imsres)//"_to_FV3_mapping."//trim(adjustl(otype))//".nc"

        print *, 'reading IMS index file', trim(IMS_ind_file) 

        inquire(file=trim(IMS_ind_file), exist=file_exists)

        if (.not. file_exists) then
          print *, 'observation_read_IMS_full error, index file does not exist', &
                 trim(IMS_ind_file) , ' exiting'
          stop 10
        endif
    
        error=nf90_open(trim(IMS_ind_file),nf90_nowrite, ncid)
        call netcdf_err(error, 'opening file: '//trim(IMS_ind_file) )
    
        allocate(IMS_index(j_ims, i_ims, 3))

        error=nf90_inq_varid(ncid, 'tile', id_var)
        call netcdf_err(error, 'error reading sncov indices id' )

        error=nf90_get_var(ncid, id_var, IMS_index(:,:,1))
        call netcdf_err(error, 'error reading sncov indices' )
    
        error=nf90_inq_varid(ncid, 'tile_i', id_var)
        call netcdf_err(error, 'error reading sncov indices id' )

        error=nf90_get_var(ncid, id_var, IMS_index(:,:,2))
        call netcdf_err(error, 'error reading sncov indices' )
    
        error=nf90_inq_varid(ncid, 'tile_j', id_var)
        call netcdf_err(error, 'error reading sncov indices id' )

        error=nf90_get_var(ncid, id_var, IMS_index(:,:,3))
        call netcdf_err(error, 'error reading sncov indices' )
  
        ! get the FV3 lat/lon
        error=nf90_inq_varid(ncid, 'lon_fv3', id_var)
        call netcdf_err(error, 'error reading lon_fv3 id' )

        error=nf90_get_var(ncid, id_var, lonFV3)
        call netcdf_err(error, 'error reading lon_fv3 data' )

        error=nf90_inq_varid(ncid, 'lat_fv3', id_var)
        call netcdf_err(error, 'error reading lat_fv3 id' )

        error=nf90_get_var(ncid, id_var, latFV3)
        call netcdf_err(error, 'error reading lat_fv3 data' )
    
        error=nf90_inq_varid(ncid, 'oro_fv3', id_var)
        call netcdf_err(error, 'error reading oro_fv3 id' )

        error=nf90_get_var(ncid, id_var, oroFV3)
        call netcdf_err(error, 'error reading oro_fv3 data' )
    
        error = nf90_close(ncid)

        ! calculate fraction of land within grid cell that is snow covered

        land_points = 0
        snow_points = 0
        scf = nodata_real

        do irow=1, j_ims
          do icol=1, i_ims

            if(IMS_flag(icol,irow) >= 0) then
              tile   = IMS_index(icol,irow,1)
              tile_i = IMS_index(icol,irow,2)
              tile_j = IMS_index(icol,irow,3)
              land_points(tile_i,tile_j,tile) = land_points(tile_i,tile_j,tile) + 1
              ! Mike - why IMS_flag here? are we expecting possible fractional input?
              snow_points(tile_i,tile_j,tile) = snow_points(tile_i,tile_j,tile) + IMS_flag(icol,irow)
            end if

          end do
        end do

        where(land_points > 0) scf = snow_points/land_points
    
        deallocate(IMS_flag)
        deallocate(IMS_index)

        return
        
 end subroutine read_IMS_onto_model_grid

!====================================
! read in the VIIRS observations and associated index file, then
! aggregate onto the model grid.

 subroutine read_VIIRS_onto_model_grid(VIIRS_obs_file, VIIRS_ind_path, viirs_threshold, &
             jdim, idim, otype, vetfcs_in, lonFV3, latFV3, oroFV3,scf)

        implicit none

        character(len=*), intent(in)   :: VIIRS_obs_file, VIIRS_ind_path
        integer, intent(in)            :: jdim, idim
        character(len=20), intent(in)  :: otype
        real, intent(in)               :: vetfcs_in(jdim,idim,6)
        real, intent(in)               :: viirs_threshold
        real, intent(out)              :: scf(jdim,idim,6)
        real, intent(out)              :: lonFV3(jdim,idim,6)
        real, intent(out)              :: latFV3(jdim,idim,6)
        real, intent(out)              :: oroFV3(jdim,idim,6)

        integer, allocatable    :: VIIRS_flag(:,:)
        real, allocatable       :: VIIRS_scf(:,:),VIIRS_scf_tmp(:,:)
        real, allocatable       :: VIIRS_ci(:,:)
        integer, allocatable    :: VIIRS_index(:,:,:)

        integer                 :: total_land_points(jdim,idim,6), valid_points(jdim,idim,6)
        real                    :: percent_valid(jdim,idim,6)
        real                    :: VIIRS_scf_sum(jdim,idim,6)
        integer                 :: vet_int

        integer                 :: error, ncid, id_var
        integer                 :: i_viirs, j_viirs, tile, tile_i, tile_j
        logical                 :: file_exists
        character(len=250)      :: VIIRS_ind_file

        integer(hsize_t), dimension(2) :: dims
        integer(hsize_t), dimension(2) :: maxdims
        character*100,parameter:: gr_name="/HDFEOS/GRIDS/VIIRS_Daily_SnowCover_CMG/Data Fields"
        character*100,parameter:: flag_field_name="Basic_QA"
        character*100,parameter:: scf_field_name="Snow_Cover"
        character*100,parameter:: ci_field_name="Clear_Index"
        integer(hid_t)         :: file_id, gr_id
        integer(hid_t)         :: scf_dspace_id, flag_dspace_id, ci_dspace_id
        integer(hid_t)         :: flag_field_id,scf_field_id,ci_field_id

        integer                :: icol, irow

        ! read VIIRS observations in
        inquire(file=trim(VIIRS_obs_file), exist=file_exists)

        if (.not. file_exists) then
           print *, 'observation_read_VIIRS_full error,file does not exist', &
                        trim(VIIRS_obs_file) , ' exiting'
           stop 10
        else
           print *, 'Starting reading VIIRS h5 file'
        endif

        ! Initialize Fortran Interface
        call h5open_f(error)

        ! Open an existing file
        call h5fopen_f(trim(VIIRS_obs_file), H5F_ACC_RDONLY_F, file_id,error)
        if (error .ne. 0) then
           print *, '[ERROR] read_VIIRS_onto_model_grid cannot open file'
           stop 100
        endif

        ! Open an existing group
        call h5gopen_f(file_id, gr_name, gr_id, error)

        call h5dopen_f(gr_id, flag_field_name, flag_field_id, error)

        call h5dopen_f(gr_id, scf_field_name, scf_field_id, error)

        call h5dopen_f(gr_id, ci_field_name, ci_field_id, error)

        ! Get space of the dataset
        call h5dget_space_f(scf_field_id, scf_dspace_id, error)

        call h5dget_space_f(flag_field_id, flag_dspace_id, error)

        call h5dget_space_f(ci_field_id, ci_dspace_id, error)

        ! Size of the arrays
        call h5sget_simple_extent_dims_f(scf_dspace_id, dims, maxdims, error)
        if (error .eq. -1) then
           print *, '[ERROR] read_VIIRS_onto_model_grid cannot get scf size'
           stop 100
        endif

        j_viirs = dims(1)
        i_viirs = dims(2)
        allocate(VIIRS_flag(j_viirs, i_viirs))
        allocate(VIIRS_ci(j_viirs, i_viirs))
        allocate(VIIRS_scf(j_viirs,i_viirs))
        allocate(VIIRS_scf_tmp(j_viirs,i_viirs))

        ! Read H5 VIIRS file
        call h5dread_f(scf_field_id, H5T_NATIVE_DOUBLE, VIIRS_scf_tmp, dims, error,&
                       H5S_ALL_F, scf_dspace_id )
        if (error .ne. 0) then
          print *, '[ERROR] read_VIIRS_onto_model_grid cannot extract scf'
          stop 100
        endif
        
         call h5dread_f(flag_field_id,H5T_NATIVE_INTEGER,VIIRS_flag, dims, error, &
                       H5S_ALL_F, flag_dspace_id )
        if (error .ne. 0) then
          print *, '[ERROR] read_VIIRS_onto_model_grid cannot extract flag'
          stop 100
        endif

        call h5dread_f(ci_field_id,H5T_NATIVE_DOUBLE,VIIRS_ci, dims, error,&
                       H5S_ALL_F, ci_dspace_id )
        if (error .ne. 0) then
          print *, '[ERROR] read_VIIRS_onto_model_grid cannot extract ci'
          stop 100
        endif

        ! Close H5 data fields
        call h5dclose_f(flag_field_id, error)

        call h5dclose_f(scf_field_id, error)

        call h5dclose_f(ci_field_id, error)

        ! Close H5 group
        call h5gclose_f(gr_id, error)

        ! Close H5 file
        call h5fclose_f(file_id, error)

        ! Close H5 fortran interface
        call h5close_f(error)
        
        ! VIIRS codes: 0 = good,
        !            : 1 = poor,
        !            : 2 = bad,
        !            : 3 = other,
        !            : 201 211 237 239 243 250 251 252 253 254
        !            : no_decision night lake ocean &
        !              Antarctica cloud missing_L1B_data cal_fail_L1B_data bowtie_trim L1B_fill'
        where(VIIRS_scf_tmp .gt. 100.) VIIRS_scf_tmp = nodata_real
        where(VIIRS_flag .gt. 1 ) VIIRS_scf_tmp = nodata_real ! retain good/poor qc/qa, set outside range to NA
        where(VIIRS_ci .lt. 50. .or. VIIRS_ci .gt. 100.) VIIRS_scf_tmp = nodata_real ! retain clearer views without too much cloud


        ! read index file for mapping VIIRS to model grid

        VIIRS_ind_file = trim(VIIRS_ind_path)//"VIIRS_to_FV3_mapping."//trim(adjustl(otype))//".nc"

        print *, 'reading VIIRS index file', trim(VIIRS_ind_file)

        inquire(file=trim(VIIRS_ind_file), exist=file_exists)

        if (.not. file_exists) then
          print *, 'observation_read_VIIRS_full error, index file does not exist', &
                 trim(VIIRS_ind_file) , ' exiting'
          stop 10
        endif

        error=nf90_open(trim(VIIRS_ind_file),nf90_nowrite, ncid)
        call netcdf_err(error, 'opening file: '//trim(VIIRS_ind_file) )

        allocate(VIIRS_index(j_viirs, i_viirs, 3))

        error=nf90_inq_varid(ncid, 'tile', id_var)
        call netcdf_err(error, 'error reading sncov indices id' )

        error=nf90_get_var(ncid, id_var, VIIRS_index(:,:,1))
        call netcdf_err(error, 'error reading sncov indices' )

        error=nf90_inq_varid(ncid, 'tile_i', id_var)
        call netcdf_err(error, 'error reading sncov indices id' )

        error=nf90_get_var(ncid, id_var, VIIRS_index(:,:,2))
        call netcdf_err(error, 'error reading sncov indices' )

        error=nf90_inq_varid(ncid, 'tile_j', id_var)
        call netcdf_err(error, 'error reading sncov indices id' )

        error=nf90_get_var(ncid, id_var, VIIRS_index(:,:,3))
        call netcdf_err(error, 'error reading sncov indices' )

        ! get the FV3 lat/lon
        error=nf90_inq_varid(ncid, 'lon_fv3', id_var)
        call netcdf_err(error, 'error reading lon_fv3 id' )

        error=nf90_get_var(ncid, id_var, lonFV3)
        call netcdf_err(error, 'error reading lon_fv3 data' )

        error=nf90_inq_varid(ncid, 'lat_fv3', id_var)
        call netcdf_err(error, 'error reading lat_fv3 id' )

        error=nf90_get_var(ncid, id_var, latFV3)
        call netcdf_err(error, 'error reading lat_fv3 data' )

        error=nf90_inq_varid(ncid, 'oro_fv3', id_var)
        call netcdf_err(error, 'error reading oro_fv3 id' )

        error=nf90_get_var(ncid, id_var, oroFV3)
        call netcdf_err(error, 'error reading oro_fv3 data' )

        error = nf90_close(ncid)

        ! flip matrix to align with mapping index
        do irow = 1, j_viirs
            do icol = 1, i_viirs
                !convert to fractional snow cover fraction
                VIIRS_scf(irow, icol) = VIIRS_scf_tmp(irow, i_viirs-icol+1)/100.
            end do
        end do

        ! Map onto fv3 grid
        total_land_points = 0
        valid_points = 0
        VIIRS_scf_sum = 0.
        percent_valid = nodata_real
        scf = nodata_real

        do irow=1, j_viirs !7200
          do icol=1,i_viirs !3600
              tile   = VIIRS_index(irow,icol,1)
              tile_i = VIIRS_index(irow,icol,2)
              tile_j = VIIRS_index(irow,icol,3)
              vet_int = int(vetfcs_in(tile_i,tile_j,tile))
              if(vet_int .gt. 0)then
                total_land_points(tile_i,tile_j,tile) = total_land_points(tile_i,tile_j,tile)+1
                if(VIIRS_scf(irow,icol).ge.0.)then
                  valid_points(tile_i,tile_j,tile) = valid_points(tile_i,tile_j,tile)+1
                  VIIRS_scf_sum(tile_i,tile_j,tile) = VIIRS_scf_sum(tile_i,tile_j,tile)+VIIRS_scf(irow,icol)
                end if
              end if
          end do
        end do

        where(total_land_points > 0) percent_valid = real(valid_points)/real(total_land_points)
        !viirs_threshold is read from nml based on user selection, tests were done for 0.1, 0.3, 0.5
        where(percent_valid > viirs_threshold) scf = VIIRS_scf_sum/valid_points

        deallocate(VIIRS_flag)
        deallocate(VIIRS_scf)
        deallocate(VIIRS_scf_tmp)
        deallocate(VIIRS_ci)
        deallocate(VIIRS_index)

        return

 end subroutine read_VIIRS_onto_model_grid

 subroutine netcdf_err( err, string )
    
    !--------------------------------------------------------------
    ! if a netcdf call returns an error, print out a message
    ! and stop processing.
    !--------------------------------------------------------------
    
        implicit none
    
        integer, intent(in) :: err
        character(len=*), intent(in) :: string
        character(len=80) :: errmsg
    
        if( err == nf90_noerr )return
        errmsg = nf90_strerror(err)
        print*,''
        print*,'fatal error: ', trim(string), ': ', trim(errmsg)
        print*,'stop.' 
        stop 10

        return
 end subroutine netcdf_err

!====================================
! calculate snow density from forecast fields.
! density = SWE/SND where snow present. 
!         = average from snow forecasts over land, where snow not present

 subroutine calc_density(idim, jdim, lsm, landmask, swe, snd, stc, density)

       implicit none 

       integer, intent(in) :: idim, jdim, lsm
       integer, intent(in) :: landmask(idim,jdim,6)
       real, intent(in)    :: swe(idim,jdim,6), snd(idim,jdim,6), stc(idim,jdim,6)
       real, intent(out)   :: density(idim,jdim,6)

       real :: dens_mean
       integer :: i,j,t

        ! density = swe/snd
        do t=1, 6
         do i =1,idim
          do j=1,jdim
                if (snd(i,j,t) > 0.01 ) then
                        density(i,j,t) = swe(i,j,t)/snd(i,j,t)
                elseif (snd(i,j,t) <= 0.01 .and. lsm==2) then ! for noah-MP, calculate from stc
                ! divide by 1000 as noah-MP has snd in m
                        density(i,j,t) = max(80.0,min(120.,67.92+51.25*exp((stc(i,j,t)-273.15)/2.59)))/1000.
                endif
          enddo 
         enddo
        enddo
        where (density < 0.0001) density = 0.08 ! CSD units here are wrong. Update this once model is updated.nn

        if (lsm==1) then ! for noah, use mean density 
            ! calculate mean density over land
            if (count (landmask==1 .and. snd> 0.01) > 0) then
                    ! mean density over snow-covered land
                    dens_mean = sum(density, mask = (landmask==1 .and. snd>0.01 )) &
                             / count (landmask==1 .and. snd> 0.01)
                    print *, 'mean density: ', dens_mean
            else
                    dens_mean = 0.1  ! default value if have no snow 
                    print *, 'no snow, using default density: ', dens_mean
            endif

            ! for grid cells with no valid density, fill in the average snodens
            where( snd <= 0.01 ) density = dens_mean
        endif

 end subroutine calc_density

!====================================
! calculate SWE from fractional snow cover, using the noah model relationship 
! uses empirical inversion of snow depletion curve in the model 

 subroutine calcSWE_noah(scf, vetfcs_in, swefcs, idim, jdim, swe)
        
        implicit none
        !
        integer, intent(in)     :: idim,jdim
        real, intent(in)        :: vetfcs_in(idim,jdim,6)
        real, intent(in)        :: swefcs(idim,jdim,6)
        real, intent(inout)     :: scf(idim,jdim,6) 
        real, intent(out)       :: swe(idim,jdim,6)
        
        integer            :: vetfcs(idim,jdim,6)
        !snup_array is the swe (mm) at which scf reaches 100% 
        real               :: snupx(30), snup, salp, rsnow
        integer            :: i,j,t,vtype_int


        ! fill background values to nan
        swe = nodata_real
   
        ! note: this is an empirical inversion of   snfrac rotuine in noah 
        !  should really have a land model check in here. 

        !this is for the igbp veg classification scheme.
        ! swe at which snow cover reaches 100%, in m
        snupx = (/0.080, 0.080, 0.080, 0.080, 0.080, 0.020,     &
                0.020, 0.060, 0.040, 0.020, 0.010, 0.020,                       &
                0.020, 0.020, 0.013, 0.013, 0.010, 0.020,                       &
                0.020, 0.020, 0.000, 0.000, 0.000, 0.000,                       &
                0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
    
        salp = -4.0
        vetfcs = nint(vetfcs_in)
        ! this is done in the noaa code, but we don't want to create a value outside land
        ! where(vetfcs==0) vetfcs = 7 
       
        do t =1, 6 
          do i = 1, idim 
            do j = 1, jdim 
          
                if ( abs( scf(i,j,t) - nodata_real ) > nodata_tol ) then  ! is have IMS/VIIRS data
                    if  (vetfcs(i,j,t)>0)  then ! if model has land
                        snup = snupx(vetfcs(i,j,t))*1000. ! convert to mm
                        if (snup == 0.) then
                            print*, " 0.0 snup value, check vegclasses", vetfcs(i,j,t)
                            stop 10
                        endif

                        ! if model and IMS/VIIRS both have 100% snow cover, don't convert IMS/VIIRS to a snow depth
                        if ( (swefcs(i,j,t) >= snup)  .and. (scf(i,j,t) >= 1.0 ) ) cycle 

                        if (scf(i,j,t) >= 1.0) then
                            rsnow = 1.
                        elseif (scf(i,j,t) < 0.001) then
                            rsnow = 0.0 
                        else
                            rsnow = min(log(1. - scf(i,j,t)) / salp, 1.0) 
                        endif  
                        ! return swe in mm 
                        swe(i,j,t) = rsnow * snup  !  mm
                    else  ! if model is not land, remove the IMS/VIIRS data
                        scf(i,j,t) = nodata_real 
                    endif
                endif
            enddo  
          enddo  
        enddo  
        return
    
 end subroutine calcSWE_noah

!====================================
! calculate IMS/VIIRS SD from fractional IMS/VIIRS snow cover, using the noah-MP model relationship
! (this is the inverse of calcSCF_noahmp). 

 subroutine calcSD_noahmp(scf, vetfcs_in, denfcs, idim, jdim,  snd)

        implicit none
        !
        integer, intent(in)     :: idim,jdim 
        real, intent(in)        :: vetfcs_in(idim,jdim,6)
        real, intent(in)        :: denfcs(idim,jdim,6)
        real, intent(inout)     :: scf(idim,jdim,6)
        real, intent(out)       :: snd(idim,jdim,6)

        integer            :: vetfcs
        real               :: mfsno, scffac,  bdsno, fmelt
        integer            :: i,j,t,vtype_int


        ! fill background values to nan
        snd = nodata_real

        do t =1, 6
          do i = 1, idim
            do j = 1, jdim
                if ( abs( scf(i,j,t) - nodata_real ) > nodata_tol ) then  ! if have IMS/VIIRS data
                    vetfcs = int(vetfcs_in(i,j,t))
                    if  (vetfcs>0)  then ! if model has land
                       if ( scf(i,j,t) < 0.5 ) then  ! if obs SCF< 0.5, reduce to 0.
                        snd(i,j,t) = 0. 
                       else
                        ! calculate snow depth
                        mfsno  =  mfsno_table(vetfcs)
                        scffac = scffac_table(vetfcs)
                        bdsno   = max(50., min(650.,denfcs(i,j,t)*1000.) )   ! x1000, as noah-mp has SND in m.
                        fmelt    = (bdsno/100.)**mfsno
                        snd(i,j,t) =  (scffac * fmelt)*atanh(trunc_scf)*1000. ! x1000 into mm 
                      endif
                    else 
                        scf(i,j,t) = nodata_real ! also remove SCF if not land in model.
                    endif
                endif
            enddo
          enddo
        enddo

 end subroutine calcSD_noahmp

!====================================
! calculate fractional snow cover from SD and density for Noah-MP
! copied from module_sf_noahmplsm.f90 

 subroutine calcSCF_noahmp(vetfcs_in, denfcs, sndfcs, idim, jdim, scffcs) 

        implicit none
        !
        integer, intent(in)     :: idim,jdim 
        real, intent(in)        :: vetfcs_in(idim,jdim,6)
        real, intent(in)        :: denfcs(idim,jdim,6)
        real, intent(in)        :: sndfcs(idim,jdim,6)
        real, intent(out)       :: scffcs(idim,jdim,6)

        integer            :: vetfcs
        real               :: mfsno, scffac,  bdsno, fmelt, snowh
        integer            :: i,j,t,vtype_int

        do t =1, 6
          do i = 1, idim
            do j = 1, jdim
                vetfcs = int(vetfcs_in(i,j,t))
                if  (vetfcs>0)  then ! if model has land
                     mfsno  =  mfsno_table(vetfcs)
                     scffac = scffac_table(vetfcs)
                     snowh  = sndfcs(i,j,t)*0.001 ! into m
                     if(sndfcs(i,j,t) .gt.0.0)  then
                         bdsno    = denfcs(i,j,t)*1000. ! 1000, as noah-mp has snd in m 
                         fmelt    = (bdsno/100.)**mfsno
                         scffcs(i,j,t) = tanh( snowh /(scffac * fmelt))
                     else 
                         scffcs(i,j,t) = 0.
                     endif 
               endif 
            enddo 
          enddo 
        enddo 


 end subroutine calcSCF_noahmp

 end module SCFaggregate_mod
 
