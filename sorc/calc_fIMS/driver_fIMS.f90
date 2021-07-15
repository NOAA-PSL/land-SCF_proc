 program driver_fIMS

! calculate fractional snow cover on model grid, from IMS snow cover obs. 
! then calculate SWE from fractional snow cover, assuming the snow 
! depletion curve used by the noah model. 
! 
! Clara Draper, July 2021 (based on code from Tseganeh Gichamo, Youlong Xia)


 use IMSaggregate_mod

 implicit none

 include 'mpif.h'

 integer :: idim, jdim, num_tiles
 integer :: lensfc, ierr
 integer :: nprocs, myrank, io

 character(len=10)   :: date_str ! yyyymmddhh
 character(len=500)  :: IMS_snowcover_path, IMS_indexes_path
 logical             :: file_exists

 namelist/fIMS_nml/  idim, jdim, date_str, IMS_snowcover_path, IMS_indexes_path

 data num_tiles/6/

 call mpi_init(ierr)
 call mpi_comm_size(mpi_comm_world, nprocs, ierr)
 call mpi_comm_rank(mpi_comm_world, myrank, ierr)

 print*,"starting fIMS program on rank ", myrank
 
 ! check 6 procs assigned.
 if (nprocs .NE. num_tiles) then 
        print *, 'driver_fIMS error: nprocs must = num_tiles', nprocs, num_tiles 
        call mpi_abort(mpi_comm_world, 999)
 endif

 ! read namelist
 inquire(file='fims.nml', exist=file_exists)

 if (.not. file_exists) then
        print *, 'namelistfile does not exist, exiting' 
        call mpi_abort(mpi_comm_world, 10)
 endif

 open (action='read', file='fims.nml', iostat=ierr, newunit=io)
 read (nml=fIMS_nml, iostat=ierr, unit=io) 
 close (io) 

 lensfc = idim*jdim ! total number of points for the cubed-sphere tile

 call mpi_barrier(mpi_comm_world, ierr)

 call calculate_IMS_fsca(num_tiles, myrank, idim, jdim, &
                      lensfc, date_str, IMS_snowcover_path, IMS_indexes_path)

 print*,"calcfIMS returned on rank", myrank

 call mpi_barrier(mpi_comm_world, ierr)
 call mpi_finalize(ierr)

 end program driver_fIMS
