module IMSaggregate_mod

use netcdf

private
public calculate_IMS_fsca

real, parameter    ::  nodata_real = -999. 
integer, parameter ::  nodata_int = -999
real, parameter    ::  nodata_tol = 0.1

contains

!====================================
! main routine to read in inputs, calculate IMS snow cover fraction, IMS SWE, 
! then IMS SND, and write out results on model grid.
! SWE is calculated using the model relationship. 
! SD is calculated using the forecast snow density. 
! SD is QC'ed out where both model and obs have 100% snow cover 
! (since can get no info from IMS snow cover in this case)

subroutine calculate_IMS_fsca(idim, jdim, date_str, IMS_snowcover_path, IMS_index_path)
                                                        
        implicit none
        
        integer, intent(in)            :: idim, jdim
        character(len=10), intent(in)  :: date_str ! yyyyddd
        character(len=*), intent(in)   :: IMS_snowcover_path, IMS_index_path

        character(len=250)  :: IMS_sncov_file
        real                :: IMS_sncov(idim,jdim,6)   ! IMS fractional snow cover in model grid

!=============================================================================================
! 1. Read forecast info, and IMS data and indexes from file, then calculate SWE
!=============================================================================================

        ! read IMS obs, and indexes, map to model grid
        IMS_sncov_file = trim(IMS_snowcover_path)//"ims"//trim(date_str)//"_4km_v1.3.asc"  

        print *, 'reading IMS snow cover data from ', trim(IMS_sncov_file) 

        call read_IMS_onto_model_grid(IMS_sncov_file, IMS_index_path, jdim, idim, IMS_sncov)

!=============================================================================================
! 2.  Write outputs
!=============================================================================================
        
        call write_fsca_outputs(idim, jdim, IMS_sncov) 

        return

 end subroutine calculate_IMS_fsca

!====================================
! routine to write the output to file

 subroutine write_fsca_outputs(idim, jdim, IMS_sncov)

      !------------------------------------------------------------------
      !------------------------------------------------------------------
      implicit none

      integer, intent(in)         :: idim, jdim
      real, intent(in)            :: IMS_sncov(jdim,idim,6)

      character(len=250)          :: output_file
      character(len=1)            :: tile_str
      integer                     :: fsize=65536, inital=0
      integer                     :: header_buffer_val = 16384
      integer                     :: dIMS_3d(3), dIMS_strt(3), dIMS_end(3)
      integer                     :: error, i, ncid
      integer                     :: dim_x, dim_y, dim_time
      integer                     :: id_x, id_y, id_time
      integer                     :: id_IMScov, id_IMSsnd 
      integer                     :: itile
 
      real(kind=4)                :: times
      real(kind=4), allocatable   :: x_data(:), y_data(:)
      real(kind=8), allocatable   :: dum2d(:,:)

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

        error = nf90_def_var(ncid, 'IMSfsca', nf90_double, dIMS_3d, id_IMScov)
        call netcdf_err(error, 'defining IMSfsca' )
        error = nf90_put_att(ncid, id_IMScov, "long_name", "IMS fractional snow covered area")
        call netcdf_err(error, 'defining IMSfsca long name' )
        error = nf90_put_att(ncid, id_IMScov, "units", "-")
        call netcdf_err(error, 'defining IMSfsca units' )

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
        
        error = nf90_put_var(ncid, id_IMScov, IMS_sncov(:,:,itile), dIMS_strt, dIMS_end)
        call netcdf_err(error, 'writing IMSfsca record')

        deallocate(x_data, y_data)
        deallocate(dum2d)

        error = nf90_close(ncid)

      end do
    
 end subroutine write_fsca_outputs

!====================================
! read in the IMS observations and  associated index file, then 
! aggregate onto the model grid.

 subroutine read_IMS_onto_model_grid(IMS_sncov_file, IMS_index_path, &
                    jdim, idim, IMS_sncov)
                    
        implicit none
    
        character(len=*), intent(in)   :: IMS_sncov_file, IMS_index_path
        integer, intent(in)            :: jdim, idim 
        real, intent(out)              :: IMS_sncov(jdim,idim,6)     
    
        integer, allocatable    :: IMS_flag(:,:)   
        integer, allocatable    :: IMS_index(:,:,:)
        real                    :: land_points(jdim,idim,6), snow_points(jdim,idim,6)
        
        integer                :: error, ncid, id_dim, id_var , n_ind
        integer                :: i_ims, j_ims, itile, tile, tile_i, tile_j
        logical                :: file_exists
        character(len=3)       :: resl_str
        character(len=250)     :: IMS_index_file

        integer                :: icol, irow

        ! read IMS observations in
        inquire(file=trim(IMS_sncov_file), exist=file_exists)

        if (.not. file_exists) then
           print *, 'observation_read_IMS_full error,file does not exist', &
                        trim(IMS_sncov_file) , ' exiting'
           stop
        endif

        i_ims = 6144
        j_ims = 6144
        allocate(IMS_flag(j_ims, i_ims))   
          
        open(10, file=IMS_sncov_file, form="formatted", status="old")
          
        do irow = 1, 30
          read(10,*)     ! read ims ascii header
        end do

        do irow = 1, j_ims
          read(10,'(6144i1)') (IMS_flag(icol, irow), icol=1, i_ims)
        end do

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

        write(resl_str, "(i3)") idim

        IMS_index_file = trim(IMS_index_path)//"fv3_mapping_C"//trim(adjustl(resl_str))//".nc"                       
        print *, 'reading IMS index file', trim(IMS_index_file) 

        inquire(file=trim(IMS_index_file), exist=file_exists)

        if (.not. file_exists) then
          print *, 'observation_read_IMS_full error, index file does not exist', &
                 trim(IMS_index_file) , ' exiting'
          stop
        endif
    
        error=nf90_open(trim(IMS_index_file),nf90_nowrite, ncid)
        call netcdf_err(error, 'opening file: '//trim(IMS_index_file) )
    
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
    
        error = nf90_close(ncid)

        ! calculate fraction of land within grid cell that is snow covered

        land_points = 0
        snow_points = 0
        IMS_sncov = nodata_real

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

        where(land_points > 0) IMS_sncov = snow_points/land_points
    
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
 
