module IMSaggregate_mod

use netcdf

private
public calculate_scfIMS

real, parameter    ::  nodata_real = -999. 
integer, parameter ::  nodata_int = -999
real, parameter    ::  nodata_tol = 0.1

! IMS noah-MP snow depth retrieval parameters
real, parameter :: trunc_scf = 0.95 ! For the Noah-MP snow depletion curve, SCF asymptotes to 1. as SD increases 
                                    ! use this value when calculating SD to represent "full" coverage
real, parameter :: sndIMS_max = 300. ! maximum snow depth for snow depth derived from SCF.


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
! main routine to read in inputs, calculate IMS snow cover fraction, IMS SWE, 
! then IMS SND, and write out results on model grid.
! SWE is calculated using the model relationship. 
! SD is calculated using the forecast snow density. 
! SD is QC'ed out where both model and obs have "full" snow cover 
! (since can get no info from IMS snow cover in this case)

! note: original version of this code used model and observed snow cover as a fraction. 
! However, comparison between model and observations shows that the obs can have <100% 
! snow cover for very deep snow (likely, due to satellite seeing trees, etc), while this 
! does not (cannot) occur in the model.  This resulted in the DA incorrectly removing snow 
! when this occured. The code has now been reverted to process the IMS SCF as 100% for SFC>50%, and 
! as 0% for SCF < 50%.

! note on timing: files are once a day, nominally at 0 UTC. 
subroutine calculate_scfIMS(idim, jdim, otype, yyyymmddhh, jdate, IMS_obs_path, & 
                            IMS_ind_path, fcst_path, lsm, imsformat, imsversion, imsres, skip_SD)
                                                        
        implicit none
        
        integer, intent(in)            :: idim, jdim, lsm
        integer, intent(in)           :: imsformat
        character(len=20), intent(in) :: otype  
        character(len=11), intent(in)  :: yyyymmddhh
        character(len=7), intent(in)  :: jdate
        character(len=*), intent(in)   :: IMS_obs_path, IMS_ind_path, fcst_path
        logical, intent(in)            :: skip_SD
        character(len=10), intent(in)  :: imsversion, imsres 

        real                :: vtype(idim,jdim,6)       ! model vegetation type
        integer             :: landmask(idim,jdim,6)
        real                :: swefcs(idim,jdim,6), sndfcs(idim,jdim,6) ! forecast SWE, SND
        real                :: stcfcs(idim,jdim,6) ! forecast soil temp  
        real                :: denfcs(idim,jdim,6) ! forecast density
        real                :: scffcs(idim,jdim,6) ! forecast snow cover fraction
        real                :: scfIMS(idim,jdim,6) ! IMS snow cover fraction, on model grid
        real                :: sweIMS(idim,jdim,6) ! SWE derived from scfIMS, on model grid - only needed for noah lsm
        real                :: sndIMS(idim,jdim,6) ! snow depth derived from scfIMS, on model grid
        real                :: lonFV3(idim,jdim,6) ! longitude on model grid
        real                :: latFV3(idim,jdim,6) ! latutide on model grid
        real                :: oroFV3(idim,jdim,6) ! orography model grid
        character(len=250)  :: IMS_obs_file
        integer             :: i,j,t

!=============================================================================================
! 1. Read forecast info, and IMS data and indexes from file, then calculate SWE
!=============================================================================================

        if (.not. skip_SD) then
            ! note: snow cover being read here is calculated at the start of the 
            !       time step, and does not account for snow changes over the last 
            !       time step. Resulting errors are generally very small.
            call  read_fcst(fcst_path, yyyymmddhh, idim, jdim, vtype, swefcs,  & 
                            sndfcs, stcfcs, landmask)

            call calc_density(idim, jdim, lsm, landmask, swefcs, sndfcs, stcfcs, denfcs)
        endif 

        ! read in either ascii or nc  IMS obs, and indexes, map to model grid

        if (imsformat==1) then
            IMS_obs_file = trim(IMS_obs_path)//"ims"//trim(jdate)//"_"//trim(imsres)//"_v"//trim(imsversion)//".asc"
        elseif (imsformat==2) then  
            IMS_obs_file = trim(IMS_obs_path)//"ims"//trim(jdate)//"_"//trim(imsres)//"_v"//trim(imsversion)//".nc"
        else
          print *, 'fatal error reading IMS snow cover data '    
       endif
       
       print *, 'reading IMS snow cover data from ', trim(IMS_obs_file) 

       call read_IMS_onto_model_grid(IMS_obs_file, IMS_ind_path, imsformat, imsres, &
                                   jdim, idim, otype, lonFV3, latFV3, oroFV3, scfIMS)
       
       if (.not. skip_SD) then
            ! calculate SWE from IMS snow cover fraction (using model relationship)
            ! no value is calculated if both IMS and model have 100% snow cover
            ! also removes scfIMS if model is non-land
            if (lsm==1) then
                call calcSWE_noah(scfIMS, vtype, swefcs, idim, jdim, sweIMS)
                ! calculate snow depth from IMS SWE, using model density
                do t=1,6
                  do i=1,idim
                    do j=1,jdim 
                      if  ( abs( sweIMS(i,j,t) -nodata_real ) > nodata_tol ) then
                        sndIMS(i,j,t) = sweIMS(i,j,t)/denfcs(i,j,t)
                      else
                        sndIMS(i,j,t) = nodata_real
                      endif
                    enddo
                  enddo
                enddo
            elseif (lsm==2) then 
                ! calculate SD from IMS SCF
                call calcSD_noahmp(scfIMS, vtype, denfcs, idim, jdim, sndIMS)

                ! calculate SCF from model forecast SD and SWE (since not always in restart)
                call calcSCF_noahmp(vtype, denfcs, sndfcs, idim, jdim, scffcs) 
                
                !exclude IMS snow depth, where both IMS and model are 100% 
                do t=1,6
                  do i=1,idim
                    do j=1,jdim 
                            ! do not assimilate, if obs and model indicate "full" snow. 
                            ! (recall: converting obs SCF>0.5 to full snow)
                            ! additional check to remove obs if obs have full snow, but calculated snow depth is lower than the model
                            if ( (scfIMS(i,j,t) >= 0.5 ) .and.  & 
                                    (  (scffcs(i,j,t) > trunc_scf ) .or.  (sndfcs(i,j,t ) >  sndIMS(i,j,t) ) ) ) then 
                                   sndIMS(i,j,t) = nodata_real
                            endif 
                            if ( sndIMS(i,j,t) > sndIMS_max ) then 
                                   sndIMS(i,j,t) = nodata_real
                            endif
                    enddo 
                  enddo 
                enddo

            else 
                print *, 'unknown lsm:', lsm, ', choose 1 - noah, 2 - noah-mp' 
                stop 10 
            endif
        else 
                sndIMS=nodata_real
        endif ! skip_SD
!=============================================================================================
! 2.  Write outputs
!=============================================================================================
        
        !call write_IMS_outputs_2D(idim, jdim, scfIMS,sndIMS)
        call write_IMS_outputs_vec(idim, jdim, otype, yyyymmddhh, scfIMS, sndIMS, lonFV3, latFV3, oroFV3)
        
        return

 end subroutine calculate_scfIMS

!====================================
! routine to write the output to file - 2D (one file per tile)

 subroutine write_IMS_outputs_2D(idim, jdim, scfIMS, sndIMS)

      !------------------------------------------------------------------
      !------------------------------------------------------------------
      implicit none

      integer, intent(in)         :: idim, jdim
      real, intent(in)            :: scfIMS(idim,jdim,6)
      real, intent(in)            :: sndIMS(idim,jdim,6)

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
        
        error = nf90_put_var(ncid, id_scfIMS, scfIMS(:,:,itile), dIMS_strt, dIMS_end)
        call netcdf_err(error, 'writing IMSscf record')

        error = nf90_put_var(ncid, id_sndIMS, sndIMS(:,:,itile), dIMS_strt, dIMS_end)
        call netcdf_err(error, 'writing IMSsnd record')

        error = nf90_close(ncid)

      end do
    
 end subroutine write_IMS_outputs_2D


!====================================
! routine to write the output to file - vector
! also writes out the model lat/lon for the grid cell that the data have been 
! processed onto.

 subroutine write_IMS_outputs_vec(idim, jdim, otype, date_str, scfIMS, sndIMS, lonFV3, latFV3, oroFV3)

    implicit none
    character(len=11), intent(in)  :: date_str
    character(len=20), intent(in)  :: otype
    integer, intent(in)         :: idim, jdim
    real, intent(in)            :: scfIMS(idim,jdim,6)
    real, intent(in)            :: sndIMS(idim,jdim,6)
    real, intent(in)            :: latFV3(idim,jdim,6)
    real, intent(in)            :: lonFV3(idim,jdim,6)
    real, intent(in)            :: oroFV3(idim,jdim,6)

    double precision            :: sec_since
    integer                     :: offset_ss
    integer                     :: dim_time, id_time

    character(len=4)            :: yyyy
    character(len=2)            :: mm, dd
    character(len=19)           :: current_date
    character(len=19)           :: since_date = "1970-01-01 00:00:00"

    character(len=250)          :: output_file
    integer                     :: header_buffer_val = 16384
    integer                     :: i,j,t,n, nobs
    integer                     :: error, ncid
    integer                     :: id_scfIMS, id_sndIMS , id_obs, id_lon, id_lat, id_oro
    real, allocatable           :: data_vec(:,:)
    real, allocatable           :: coor_vec(:,:)

    yyyy=date_str(1:4)
    mm=date_str(5:6)
    dd=date_str(7:8)

    offset_ss=0
    current_date=yyyy//'-'//mm//'-'//dd//' 18:00:00'

    call calc_sec_since(since_date, current_date, offset_ss, sec_since)
    
    output_file = "./IMSscf."//date_str(1:8)//"."//trim(adjustl(otype))//".nc"
    print*,'writing output to ',trim(output_file) 
    
    !--- create the file
    !error = nf90_create(output_file, ior(nf90_netcdf4,nf90_classic_model), ncid, initialsize=inital, chunksize=fsize)
    error = nf90_create(output_file, ior(nf90_netcdf4,nf90_classic_model), ncid)
    call netcdf_err(error, 'creating file='//trim(output_file) )

    ! collect obs
    ! note: writing out all scf ovs. snd will be missing for many locations 
    !       since it is removed if both model and ovs have scf==1
    nobs = count (abs(scfIMS -nodata_real) > nodata_tol)
    print *, 'writing out', nobs, ' observations'

    allocate(data_vec(2,nobs)) 
    allocate(coor_vec(3,nobs)) 

    !--- define location dimension
    error = nf90_def_dim(ncid, 'numobs', nobs, id_obs)
    call netcdf_err(error, 'defining obs dimension' )

    ! --- define time dimension
    error = nf90_def_dim(ncid, 'time', nf90_unlimited, dim_time)
    call netcdf_err(error, 'defining time dimension' )
 
    !--- define time
    error = nf90_def_var(ncid, 'time', nf90_double, dim_time, id_time)
    call netcdf_err(error, 'defining time' )
    error = nf90_put_att(ncid, id_time, "long_name", "time")
    call netcdf_err(error, 'defining time long name' )
    error = nf90_put_att(ncid, id_time, "units", "seconds since "//since_date)
    call netcdf_err(error, 'defining time units' ) 

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
    error = nf90_def_var(ncid, 'IMSscf', nf90_double, (/id_obs,dim_time/), id_scfIMS)
    call netcdf_err(error, 'defining IMSscf' )
    error = nf90_put_att(ncid, id_scfIMS, "long_name", "IMS snow covered fraction")
    call netcdf_err(error, 'defining IMSscf long name' )
    error = nf90_put_att(ncid, id_scfIMS, "units", "-")
    call netcdf_err(error, 'defining IMSscf units' )

    !--- define snow depth
    error = nf90_def_var(ncid, 'IMSsnd', nf90_double, (/id_obs,dim_time/), id_sndIMS)
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
        if (abs(scfIMS(i,j,t) -nodata_real) > nodata_tol) then 
                n=n+1
                data_vec(1,n) = scfIMS(i,j,t)
                data_vec(2,n) = sndIMS(i,j,t)
                coor_vec(1,n) = lonFV3(i,j,t)
                coor_vec(2,n) = latFV3(i,j,t)
                coor_vec(3,n) = oroFV3(i,j,t)
        endif
      enddo 
     enddo 
    enddo

    ! --- put time, lat, lon, data
    error = nf90_put_var(ncid, id_time, sec_since)
    call netcdf_err(error, 'writing time record')

    error = nf90_put_var(ncid, id_scfIMS, data_vec(1,:), start = (/1,1/), count = (/nobs,1/))
    call netcdf_err(error, 'writing IMSscf record')

    error = nf90_put_var(ncid, id_sndIMS, data_vec(2,:), start = (/1,1/), count = (/nobs,1/))
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
! read in required forecast fields from a UFS surface restart 

 subroutine read_fcst(path, date_str, idim, jdim, vetfcs, swefcs, sndfcs, stcfcs, landmask)

        implicit none
        character(len=*), intent(in)      :: path
        character(11), intent(in)          :: date_str
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
            fcst_file = trim(path)//trim(date_str)//"0000.sfc_data.tile"//tt//".nc"
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
! read in the IMS observations and  associated index file, then 
! aggregate onto the model grid.

 subroutine read_IMS_onto_model_grid(IMS_obs_file, IMS_ind_path, &
            imsformat, imsres, jdim, idim, otype, lonFV3, latFV3, oroFV3,scfIMS)
                    
        implicit none
    
        character(len=*), intent(in)   :: IMS_obs_file, IMS_ind_path
        integer, intent(in)            :: jdim, idim, imsformat
        character(len=20), intent(in)  :: otype
        character(len=10), intent(in)  :: imsres
        real, intent(out)              :: scfIMS(jdim,idim,6)     
        real, intent(out)              :: lonFV3(jdim,idim,6)     
        real, intent(out)              :: latFV3(jdim,idim,6)     
        real, intent(out)              :: oroFV3(jdim,idim,6)     
    
        integer, allocatable    :: IMS_flag(:,:)   
        integer, allocatable    :: IMS_index(:,:,:)
        real                    :: land_points(jdim,idim,6), snow_points(jdim,idim,6)
        
        integer                :: error, ncid, id_dim, id_var , n_ind
        integer                :: i_ims, j_ims, itile, tile, tile_i, tile_j
        logical                :: file_exists
        character(len=250)     :: IMS_ind_file

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

        elseif (imsformat==2) then
        ! read in netCDF IMS data
           error=nf90_open(trim(IMS_obs_file),nf90_nowrite, ncid)
           call netcdf_err(error, 'opening file: '//trim(IMS_obs_file) )

           error=nf90_inq_varid(ncid, 'IMS_Surface_Values', id_var)
           call netcdf_err(error, 'error reading IMS id' )
           
           error=nf90_get_var(ncid, id_var, IMS_flag)
           call netcdf_err(error, 'error reading IMS nc data' )
           
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
        scfIMS = nodata_real

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
        
        where(land_points > 0) scfIMS = snow_points/land_points
    
        deallocate(IMS_flag)
        deallocate(IMS_index)
        
        return
        
 end subroutine read_IMS_onto_model_grid

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

 subroutine calcSWE_noah(scfIMS, vetfcs_in, swefcs, idim, jdim, sweIMS)
        
        implicit none
        !
        integer, intent(in)     :: idim,jdim
        real, intent(in)        :: vetfcs_in(idim,jdim,6)
        real, intent(in)        :: swefcs(idim,jdim,6)
        real, intent(inout)     :: scfIMS(idim,jdim,6) 
        real, intent(out)       :: sweIMS(idim,jdim,6)
        
        integer            :: vetfcs(idim,jdim,6)
        !snup_array is the swe (mm) at which scf reaches 100% 
        real               :: snupx(30), snup, salp, rsnow
        integer            :: i,j,t,vtype_int


        ! fill background values to nan
        sweIMS = nodata_real
   
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
          
                if ( abs( scfIMS(i,j,t) - nodata_real ) > nodata_tol ) then  ! is have IMS data
                    if  (vetfcs(i,j,t)>0)  then ! if model has land
                        snup = snupx(vetfcs(i,j,t))*1000. ! convert to mm
                        if (snup == 0.) then
                            print*, " 0.0 snup value, check vegclasses", vetfcs(i,j,t)
                            stop 10
                        endif

                        ! if model and IMS both have 100% snow cover, don't convert IMS to a snow depth
                        if ( (swefcs(i,j,t) >= snup)  .and. (scfIMS(i,j,t) >= 1.0 ) ) cycle 

                        if (scfIMS(i,j,t) >= 1.0) then
                            rsnow = 1.
                        elseif (scfIMS(i,j,t) < 0.001) then
                            rsnow = 0.0 
                        else
                            rsnow = min(log(1. - scfIMS(i,j,t)) / salp, 1.0) 
                        endif  
                        ! return swe in mm 
                        sweIMS(i,j,t) = rsnow * snup  !  mm
                    else  ! if model is not land, remove the IMS data
                        scfIMS(i,j,t) = nodata_real 
                    endif
                endif
            enddo  
          enddo  
        enddo  
        return
    
 end subroutine calcSWE_noah

!====================================
! calculate IMS SD from fractional IMS snow cover, using the noah-MP model relationship
! (this is the inverse of calcSCF_noahmp). 

 subroutine calcSD_noahmp(scfIMS, vetfcs_in, denfcs, idim, jdim,  sndIMS)

        implicit none
        !
        integer, intent(in)     :: idim,jdim 
        real, intent(in)        :: vetfcs_in(idim,jdim,6)
        real, intent(in)        :: denfcs(idim,jdim,6)
        real, intent(inout)     :: scfIMS(idim,jdim,6)
        real, intent(out)       :: sndIMS(idim,jdim,6)

        integer            :: vetfcs
        real               :: mfsno, scffac,  bdsno, fmelt
        integer            :: i,j,t,vtype_int


        ! fill background values to nan
        sndIMS = nodata_real

        do t =1, 6
          do i = 1, idim
            do j = 1, jdim
                if ( abs( scfIMS(i,j,t) - nodata_real ) > nodata_tol ) then  ! if have IMS data
                    vetfcs = int(vetfcs_in(i,j,t))
                    if  (vetfcs>0)  then ! if model has land
                       if ( scfIMS(i,j,t) < 0.5 ) then  ! if obs SCF< 0.5, reduce to 0.
                        sndIMS(i,j,t) = 0. 
                       else
                        ! calculate snow depth
                        mfsno  =  mfsno_table(vetfcs)
                        scffac = scffac_table(vetfcs)
                        bdsno   = max(50., min(650.,denfcs(i,j,t)*1000.) )   ! x1000, as noah-mp has SND in m.
                        fmelt    = (bdsno/100.)**mfsno
                        sndIMS(i,j,t) =  (scffac * fmelt)*atanh(trunc_scf)*1000. ! x1000 into mm 
                      endif
                    else 
                        scfIMS(i,j,t) = nodata_real ! also remove SCF if not land in model.
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

 end module IMSaggregate_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate time in seconds
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_sec_since(since_date, current_date, offset_ss, sec_since)

! calculate number of seconds between since_date and current_date

double precision      :: sec_since
character*19 :: since_date, current_date  ! format: yyyy-mm-dd hh:nn:ss
integer      :: offset_ss
integer      :: since_yyyy, since_mm, since_dd, since_hh, since_nn, since_ss
integer      :: current_yyyy, current_mm, current_dd, current_hh, current_nn, current_ss
logical      :: leap_year = .false.
integer      :: iyyyy, imm
integer, dimension(12), parameter :: days_in_mm = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  write(*,*) "since_date=", since_date
  write(*,*) "current_date=", current_date
  write(*,*) "offset_ss=", offset_ss
  sec_since = 0

  read(since_date( 1: 4),  '(i4)') since_yyyy
  read(since_date( 6: 7),  '(i2)') since_mm
  read(since_date( 9:10),  '(i2)') since_dd
  read(since_date(12:13),  '(i2)') since_hh
  read(since_date(15:16),  '(i2)') since_nn
  read(since_date(18:19),  '(i2)') since_ss

  read(current_date( 1: 4),  '(i4)') current_yyyy
  read(current_date( 6: 7),  '(i2)') current_mm
  read(current_date( 9:10),  '(i2)') current_dd
  read(current_date(12:13),  '(i2)') current_hh
  read(current_date(15:16),  '(i2)') current_nn
  read(current_date(18:19),  '(i2)') current_ss

! not worrying about the complexity of non-recent leap years
! calculate number of seconds in all years
  do iyyyy = since_yyyy, current_yyyy
    if(mod(iyyyy,4) == 0) then
      sec_since = sec_since + 366*86400
    else
      sec_since = sec_since + 365*86400
    end if
  end do

! remove seconds from since_year
  if(mod(since_yyyy,4) == 0) leap_year = .true.

  do imm = 1,since_mm-1
    sec_since = sec_since - days_in_mm(imm)*86400
  end do

  if(leap_year .and. since_mm > 2) sec_since = sec_since - 86400

  sec_since = sec_since - (since_dd - 1) * 86400

  sec_since = sec_since - (since_hh) * 3600

  sec_since = sec_since - (since_nn) * 60

  sec_since = sec_since - (since_ss)

! remove seconds in current_year

  leap_year = .false.
  if(mod(current_yyyy,4) == 0) leap_year = .true.

  do imm = current_mm+1, 12
    sec_since = sec_since - days_in_mm(imm)*86400
  end do
  if(leap_year .and. current_mm < 3) sec_since = sec_since - 86400

  sec_since = sec_since - (days_in_mm(current_mm) - current_dd) * 86400

  sec_since = sec_since - (23 - current_hh) * 3600

  sec_since = sec_since - (59 - current_nn) * 60

  sec_since = sec_since - (60 - current_ss)

  sec_since = sec_since + offset_ss

end subroutine calc_sec_since
 
