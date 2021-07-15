module IMSaggregate_mod

use netcdf

private
public calculate_IMS_fsca

real, parameter    ::  nodata_real = -999. 
integer, parameter ::  nodata_int = -999
real, parameter    ::  nodata_tol = 0.1

contains

!====================================
! main routine to read in inputs, calculate IMS snow cover fraction, and IMS SWE, 
! write out results on model grid.

subroutine calculate_IMS_fsca(num_tiles, myrank, idim, jdim, lensfc, & 
                                date_str,IMS_snowcover_path, IMS_indexes_path)
                                                        
        !----------------------------------------------------------------------
        ! input arguments: 
        ! idim * jdim = lensfc: number of grid cells in tile = xdim * ydim   
        ! myrank: rank/id of the mpi process
        !----------------------------------------------------------------------
        implicit none
        !
        include 'mpif.h'
        
        integer, intent(in)    :: num_tiles, myrank, idim, jdim, lensfc
        character(len=10), intent(in) :: date_str ! yyyymmddhh
        character(len=*), intent(in)   :: IMS_snowcover_path, IMS_indexes_path

        integer             :: ierr     
        character(len=250)  :: IMS_inp_file, IMS_inp_file_indices 
        character(len=250)  :: IMS_out_file, vegt_inp_file
        character(len=1)    :: tile_str
        character(len=3)    :: resl_str
        real                :: sncov_IMS(lensfc)  ! IMS fractional snow cover in model grid
        real                :: swe_IMS(lensfc)    ! SWE derived from sncov_IMS, on model grid
        real                :: vetfcs(lensfc)     ! model vegetation type

        integer, parameter :: printrank = 4 


!=============================================================================================
! 1. Read veg type, and IMS data and indexes from file, then calculate SWE
!=============================================================================================

        ! rank 0 reads tile 1, etc.
        write(tile_str, '(i1.1)') (myrank+1) ! assuming <10 tiles.
        write(resl_str, "(i3)") idim

        ! read vegetation type
        vegt_inp_file = './C'//trim(adjustl(resl_str))//'.vegetation_type.tile'//tile_str//'.nc'
        if (myrank==printrank) print *, 'reading model backgroundfile for veg type', trim(vegt_inp_file) 
                                     
        call read_vegtype(vegt_inp_file, lensfc, vetfcs)

        ! read IMS obs, and indexes, map to model grid
        IMS_inp_file = trim(IMS_snowcover_path)//"IMS.SNCOV."//date_str//".nc"                 
        if (myrank==printrank) print *, 'reading IMS snow cover data from ', trim(IMS_inp_file) 

        IMS_inp_file_indices = trim(IMS_indexes_path)//"C"//trim(adjustl(resl_str))// &
                                                    ".IMS.Indices.tile"//tile_str//".nc"                       
        if (myrank==printrank) print *, 'reading IMS index file', trim(IMS_inp_file_indices) 

        call read_IMS_onto_model_grid(IMS_inp_file, IMS_inp_file_indices, &
                                                    myrank, jdim, idim, sncov_IMS)

        if (myrank==printrank) print*,'read in sncov, converting to snow depth' 
 
        call calcSWEfromSCFnoah(sncov_IMS, vetfcs, lensfc, swe_IMS)

!=============================================================================================
! 2.  Write outputs
!=============================================================================================
        
        IMS_out_file = "./IMSfSCA.tile"//tile_str//".nc"  !  
        if (myrank==printrank) print*,'writing output to ',trim(IMS_out_file) 
        call write_fsca_outputs(trim(IMS_out_file), idim, jdim, lensfc, myrank, &
                              sncov_IMS, swe_IMS) 

        call mpi_barrier(mpi_comm_world, ierr)

999 continue 
        return

 end subroutine calculate_IMS_fsca

!====================================
! routine to write the output to file

 subroutine write_fsca_outputs(output_file, idim, jdim, lensfc, myrank,   &
                                 sncov_IMS, swe_IMS)  !, anl_fsca) !updated snocov
        !------------------------------------------------------------------
        !------------------------------------------------------------------
        implicit none

        character(len=*), intent(in)      :: output_file
        integer, intent(in)         :: idim, jdim, lensfc, myrank
        real, intent(in)            ::  sncov_IMS(lensfc), swe_IMS(lensfc)

        integer                     :: fsize=65536, inital=0
        integer                     :: header_buffer_val = 16384
        integer                     :: dIMS_3d(3), dIMS_strt(3), dIMS_end(3)
        integer                     :: error, i, ncid
        integer                     :: dim_x, dim_y, dim_time
        integer                     :: id_x, id_y, id_time
        integer       :: id_IMScov, id_IMSsnd 
 
        real(kind=4)                :: times
        real(kind=4), allocatable   :: x_data(:), y_data(:)
        real(kind=8), allocatable   :: dum2d(:,:)

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

        error = nf90_def_var(ncid, 'IMSfsca', nf90_double, dIMS_3d, id_IMScov)
        call netcdf_err(error, 'defining IMSfsca' )
        error = nf90_put_att(ncid, id_IMScov, "long_name", "IMS fractional snow covered area")
        call netcdf_err(error, 'defining IMSfsca long name' )
        error = nf90_put_att(ncid, id_IMScov, "units", "-")
        call netcdf_err(error, 'defining IMSfsca units' )

        error = nf90_def_var(ncid, 'IMSswe', nf90_double, dIMS_3d, id_IMSsnd)
        call netcdf_err(error, 'defining IMSswe' )
        error = nf90_put_att(ncid, id_IMSsnd, "long_name", "IMS snow water equivalent")
        call netcdf_err(error, 'defining IMSswe long name' )
        error = nf90_put_att(ncid, id_IMSsnd, "units", "mm")
        call netcdf_err(error, 'defining IMSswe units' )

        error = nf90_enddef(ncid, header_buffer_val,4,0,4)
        call netcdf_err(error, 'defining header' )

        allocate(x_data(idim))
        do i = 1, idim
        x_data(i) = float(i)
        enddo
        allocate(y_data(jdim))
        do i = 1, jdim
        y_data(i) = float(i)
        enddo
        times = 1.0

        error = nf90_put_var( ncid, id_x, x_data)
        call netcdf_err(error, 'writing xaxis record' )
        error = nf90_put_var( ncid, id_y, y_data)
        call netcdf_err(error, 'writing yaxis record' )
        error = nf90_put_var( ncid, id_time, times)
        call netcdf_err(error, 'writing time record' )

        allocate(dum2d(idim,jdim))
        dIMS_strt(1:3) = 1
        dIMS_end(1) = idim
        dIMS_end(2) = jdim
        dIMS_end(3) = 1
        
        dum2d = reshape(sncov_IMS, (/idim, jdim/))
        error = nf90_put_var(ncid, id_IMScov, dum2d, dIMS_strt, dIMS_end)
        call netcdf_err(error, 'writing IMSfsca record')

        dum2d = reshape(swe_IMS, (/idim, jdim/))        
        error = nf90_put_var(ncid, id_IMSsnd, dum2d, dIMS_strt, dIMS_end)
        call netcdf_err(error, 'writing IMSswe record')


        deallocate(x_data, y_data)
        deallocate(dum2d)

        error = nf90_close(ncid)
    
 end subroutine write_fsca_outputs

!====================================
! read in vegetation file from a UFS surface restart 

 subroutine read_vegtype(vegt_inp_path, lensfc, vetfcs) 
    
        implicit none

        include "mpif.h"
        
        character(len=*), intent(in)      :: vegt_inp_path
        integer, intent(in)               :: lensfc
        real, intent(out)                 :: vetfcs(lensfc) 

        integer                   :: error, ncid
        integer                   :: idim, jdim, id_dim
        integer                   :: id_var

        real(kind=8), allocatable :: dummy(:,:)
        logical                   :: file_exists

        inquire(file=trim(vegt_inp_path), exist=file_exists)

        if (.not. file_exists) then 
                print *, 'read_vegtype error,file does not exist', &   
                        trim(vegt_inp_path) , ' exiting'
                call mpi_abort(mpi_comm_world, 10)
        endif
                

        error=nf90_open(trim(vegt_inp_path), nf90_nowrite,ncid)
        call netcdf_err(error, 'opening file: '//trim(vegt_inp_path) )

        error=nf90_inq_dimid(ncid, 'nx', id_dim)
        call netcdf_err(error, 'reading nx' )
        error=nf90_inquire_dimension(ncid,id_dim,len=idim)
        call netcdf_err(error, 'reading nx' )

        error=nf90_inq_dimid(ncid, 'ny', id_dim)
        call netcdf_err(error, 'reading ny' )
        error=nf90_inquire_dimension(ncid,id_dim,len=jdim)
        call netcdf_err(error, 'reading ny' )

        if ((idim*jdim) /= lensfc) then
        print*,'fatal error reading veg type: dimensions wrong.'
        call mpi_abort(mpi_comm_world, 88)
        endif

        allocate(dummy(idim,jdim))

        error=nf90_inq_varid(ncid, "vegetation_type", id_var)
        call netcdf_err(error, 'reading vegetation type id' )
        error=nf90_get_var(ncid, id_var, dummy)
        call netcdf_err(error, 'reading vegetation type' )
        vetfcs = reshape(dummy, (/lensfc/))    

        deallocate(dummy)

        error = nf90_close(ncid)
    
 end subroutine read_vegtype

!====================================
! read in the IMS observations and  associated index file, then 
! aggregate onto the model grid.

 subroutine read_IMS_onto_model_grid(inp_file, inp_file_indices, &
                    myrank, n_lat, n_lon, sncov_IMS)
                    
        implicit none
    
        include 'mpif.h'                  
    
        character(len=*), intent(in)   :: inp_file, inp_file_indices 
        integer, intent(in)            :: myrank, n_lat, n_lon 
        real, intent(out)       :: sncov_IMS(n_lat * n_lon)     
    
        integer, allocatable    :: sncov_IMS_2d_full(:,:)   
        integer, allocatable    :: data_grid_IMS_ind(:,:,:)
        real                    :: grid_dat(n_lon, n_lat)
        
        integer                :: error, ncid, id_dim, id_var , n_ind
        integer                :: dim_len_lat, dim_len_lon
        logical                :: file_exists

        ! read IMS observations in
        inquire(file=trim(inp_file), exist=file_exists)

        if (.not. file_exists) then
                print *, 'observation_read_IMS_full error,file does not exist', &
                        trim(inp_file) , ' exiting'
                call mpi_abort(mpi_comm_world, 10)
        endif
    
        error=nf90_open(trim(inp_file),nf90_nowrite,ncid)
        call netcdf_err(error, 'opening file: '//trim(inp_file) )
    
        error=nf90_inq_dimid(ncid, 'lat', id_dim)
        call netcdf_err(error, 'error reading dimension lat' )
        
        error=nf90_inquire_dimension(ncid,id_dim,len=dim_len_lat)
        call netcdf_err(error, 'error reading size of dimension lat' )
        
        error=nf90_inq_dimid(ncid, 'lon', id_dim)
        call netcdf_err(error, 'error reading dimension lon' )
    
        error=nf90_inquire_dimension(ncid,id_dim,len=dim_len_lon)
        call netcdf_err(error, 'error reading size of dimension lon' )
    

        allocate(sncov_IMS_2d_full(dim_len_lon, dim_len_lat))   

        error=nf90_inq_varid(ncid, 'Band1', id_var)
        call netcdf_err(error, 'error reading sncov id' )
        error=nf90_get_var(ncid, id_var, sncov_IMS_2d_full, start = (/ 1, 1 /), &
                                count = (/ dim_len_lon, dim_len_lat/))
        call netcdf_err(error, 'error reading sncov record' )
        
        error = nf90_close(ncid)

        ! read index file for mapping IMS to model grid 

        inquire(file=trim(inp_file_indices), exist=file_exists)

        if (.not. file_exists) then
                print *, 'observation_read_IMS_full error, index file does not exist', &
                        trim(inp_file_indices) , ' exiting'
                call mpi_abort(mpi_comm_world, 10) 
        endif
    
        error=nf90_open(trim(inp_file_indices),nf90_nowrite, ncid)
        call netcdf_err(error, 'opening file: '//trim(inp_file_indices) )
    
        error=nf90_inq_dimid(ncid, 'Indices', id_dim)
        call netcdf_err(error, 'error reading dimension indices' )
    
        error=nf90_inquire_dimension(ncid,id_dim,len=n_ind)
        call netcdf_err(error, 'error reading size of dimension indices' )

        allocate(data_grid_IMS_ind(n_ind,n_lon, n_lat))

        error=nf90_inq_varid(ncid, 'IMS_Indices', id_var)
        call netcdf_err(error, 'error reading sncov indices id' )

        error=nf90_get_var(ncid, id_var, data_grid_IMS_ind, start = (/ 1, 1, 1 /), &
                                count = (/ n_ind, n_lon, n_lat/))
        call netcdf_err(error, 'error reading sncov indices' )
    
        error = nf90_close(ncid)
   
        ! IMS codes: 0 - outside range,
        !          : 1 - sea
        !          : 2 - land, no snow 
        !          : 3 - sea ice 
        !          : 4 - snow covered land


        ! calculate fraction of land within grid cell that is snow covered
        sncov_IMS_2d_full = sncov_IMS_2d_full 
        where(sncov_IMS_2d_full == 0 ) sncov_IMS_2d_full = nodata_int ! set outside range to NA
        where(sncov_IMS_2d_full == 1 ) sncov_IMS_2d_full = nodata_int ! set sea to NA
        where(sncov_IMS_2d_full == 3 ) sncov_IMS_2d_full = nodata_int ! set sea ice to NA
        where(sncov_IMS_2d_full == 2 ) sncov_IMS_2d_full = 0          ! set land, no snow to 0
        where(sncov_IMS_2d_full == 4 ) sncov_IMS_2d_full = 1 ! set snow on land to 1

        !where(sncov_IMS_2d_full /= 4) sncov_IMS_2d_full = 0
        !where(sncov_IMS_2d_full == 4) sncov_IMS_2d_full = 1
        
        call resample_to_model_tiles_intrp(sncov_IMS_2d_full, data_grid_IMS_ind, &
                                           dim_len_lat, dim_len_lon, n_lat, n_lon, n_ind, &
                                           grid_dat)
    
        sncov_IMS = reshape(grid_dat, (/n_lat * n_lon/))
    
        deallocate(sncov_IMS_2d_full)
        deallocate(data_grid_IMS_ind)

        return
        
 end subroutine read_IMS_onto_model_grid

!====================================
! calculate SWE from fractional snow cover, using the noah model relationship 
! uses empirical inversion of snow depletion curve in the model 

 subroutine calcSWEfromSCFnoah(sncov_IMS, vetfcs_in, lensfc, swe_IMS_at_grid)
        
        implicit none
        !
        real, intent(in)        :: sncov_IMS(lensfc),  vetfcs_in(lensfc)
        integer, intent(in)     :: lensfc
        !snup_array is the swe (mm) at which scf reaches 100% 
        real, intent(out)       :: swe_IMS_at_grid(lensfc)
        
        integer            :: vetfcs(lensfc)
        real               :: snupx(30), snup, salp, rsnow
        integer                    :: indx, vtype_int


        !fill background values to nan (to differentiate those that don't have value)
        swe_IMS_at_grid = nodata_real
   
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
        !where(vetfcs==0) vetfcs = 7 
        
        do indx = 1, lensfc  
            if ( (abs(sncov_IMS(indx) -nodata_real) > nodata_tol ) .and. & 
                 (vetfcs(indx)>0)  ) then
                snup = snupx(vetfcs(indx))
                if (snup == 0.) then
                    print*, " 0.0 snup value, check vegclasses", vetfcs(indx)
                    stop
                endif

                if (sncov_IMS(indx) >= 1.0) then
                    rsnow = 1.
                elseif (sncov_IMS(indx) < 0.001) then
                    rsnow = 0.0 
                else
                    rsnow = min(log(1. - sncov_IMS(indx)) / salp, 1.0) 
                endif  
                ! return swe in mm 
                swe_IMS_at_grid(indx) = rsnow * snup * 1000. !  mm
            endif
        end do  
        return
    
 end subroutine calcSWEfromSCFnoah

!====================================
! project IMS input onto model grid
! output (grid_dat), is the average of the input (data_grid_IMS) located within each grid cell.

 subroutine resample_to_model_tiles_intrp(data_grid_IMS, data_grid_IMS_ind, &
                                            nlat_IMS, nlon_IMS, n_lat, n_lon, num_sub, & 
                                            grid_dat)
                                            
    implicit none

    integer, intent(in)     :: nlat_IMS, nlon_IMS, n_lat, n_lon, num_sub 
    integer, intent(in)     :: data_grid_IMS(nlon_IMS, nlat_IMS), data_grid_IMS_ind(num_sub, n_lon, n_lat) 
    real, intent(out)       :: grid_dat(n_lon, n_lat)

    integer   :: jc, jy, ix, num_loc_counter
    real      :: num_data
    integer   :: lonlatcoord_IMS, loncoord_IMS, latcoord_IMS
    
    grid_dat = nodata_real

    do jy=1, n_lat

        do ix=1, n_lon
            
            num_loc_counter = data_grid_IMS_ind(1, ix, jy) ! number IMS locations in this grid cell
            if (num_loc_counter < 1) then 
                cycle
            end if

            grid_dat(ix, jy) = 0.
            num_data = 0.
        
            do jc = 2, num_loc_counter+1
                lonlatcoord_IMS = data_grid_IMS_ind(jc, ix, jy) - 1 
                latcoord_IMS = lonlatcoord_IMS / nlon_IMS + 1
                loncoord_IMS = mod(lonlatcoord_IMS, nlon_IMS) + 1
                if(latcoord_IMS > nlat_IMS) then
                    latcoord_IMS = nlat_IMS
                    print*, "warning! lat coordinate outside domain boundary"
                endif
                if(loncoord_IMS > nlon_IMS) then
                    loncoord_IMS = nlon_IMS
                    print*, "warning! lon coordinate outside domain boundary"
                endif
                if ( ( data_grid_IMS(loncoord_IMS, latcoord_IMS)  - nodata_int) .ne. 0 ) then
                        grid_dat(ix, jy) =  grid_dat(ix, jy) + data_grid_IMS(loncoord_IMS, latcoord_IMS)
                        num_data = num_data+1. 
                endif
                
            end do
            if (num_data > 0 ) then
                grid_dat(ix, jy) =  grid_dat(ix, jy) / num_data
            else 
                grid_dat(ix, jy) = nodata_real
            endif

        end do

    end do

    return 

 end subroutine resample_to_model_tiles_intrp
    

 subroutine netcdf_err( err, string )
    
    !--------------------------------------------------------------
    ! if a netcdf call returns an error, print out a message
    ! and stop processing.
    !--------------------------------------------------------------
    
        implicit none
    
        include 'mpif.h'
    
        integer, intent(in) :: err
        character(len=*), intent(in) :: string
        character(len=80) :: errmsg
    
        if( err == nf90_noerr )return
        errmsg = nf90_strerror(err)
        print*,''
        print*,'fatal error: ', trim(string), ': ', trim(errmsg)
        print*,'stop.'
        call mpi_abort(mpi_comm_world, 999)
    
        return
 end subroutine netcdf_err

 end module IMSaggregate_mod
 
