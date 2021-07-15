! to call the functions here from Python, compile this file with 
! compiled with  python -m numpy.f2py -c -m imsgridmap imsgridmap.f90
! Alternatively (not working for me any more): C:/Anaconda3/Scripts/f2py.exe -c  --fcompiler=gnu95 --compiler=mingw32  -m imsgridmap imsgridmap.f90
! Then copy the .pyd and .dll to conda/envs/env_full/DLLs

                            
subroutine select_imscells_within_tile_grid(max_dx, lat_grid_ims, lon_grid_ims, &
                                   lat_grid, lon_grid, n_lat_out, n_lon_out, num_sub, &
                                   n_lat, n_lon, ntot_ims, &  !myrank, &
                                   grid_dat)
! selects IMS grid cells within a given model (FV3) grid
! Inputs:
! lat_grid_ims, lon_grid_ims: lat/lon arrays of each point on IMS file, 
! ntot_ims: size of lat_grid_ims =number of IMS grid cells= 2176*8703 for 4Km IMS
! lat_grid, lon_grid:         lat/lon arrays of coordinates of the 'corners' of the target 'grid_file'
! n_lat, n_lon = size of lat_grid, lon_grid(2*res + 1, 2*res + 1) where res = ndim in Cndim
! grid_dat = array with indices
! n_lat_out, n_lon_out, num_sub: size of grid_dat
! num_sub = maximum number of IMS grid cells within one model grid  
    Use, Intrinsic :: IEEE_ARITHMETIC
    
    Implicit None
    
    Real, Intent(in)        :: max_dx
    Integer, Intent(In)     :: ntot_ims, n_lat, n_lon, num_sub, n_lat_out, n_lon_out    !n_lon, n_lat,     !, myrank 
    Real, Intent(In)        :: lat_grid_ims(ntot_ims), lon_grid_ims(ntot_ims) 
    Real, Intent(In)        :: lat_grid(n_lat, n_lon), lon_grid(n_lat, n_lon)
    
    Integer, Intent(Out)    :: grid_dat(n_lat_out, n_lon_out, num_sub)
    
    Integer   :: j, iy, jx, num_loc_counter
    Real      :: min_y, max_y, min_x, max_x

    ! #lon_grid[lon_grid > 180.] = lon_grid[lon_grid > 180.] - 360.
    print*, "total number ims grid cells = ", ntot_ims

    Do iy=2, n_lat-1, 2
        
        !print*, "process: ", myrank, "loop ", indx

        Do jx=2, n_lon-1, 2
            print*, "iy ", iy/2, " jx ", jx/2

            min_y = MIN(lat_grid(iy-1, jx-1), lat_grid(iy-1, jx), lat_grid(iy-1, jx+1), &
                        lat_grid(iy, jx-1), lat_grid(iy, jx), lat_grid(iy, jx+1), &
                        lat_grid(iy+1, jx-1), lat_grid(iy+1, jx), lat_grid(iy+1, jx+1))
            max_y = MAX(lat_grid(iy-1, jx-1), lat_grid(iy-1, jx), lat_grid(iy-1, jx+1), &
                        lat_grid(iy, jx-1), lat_grid(iy, jx), lat_grid(iy, jx+1), &
                        lat_grid(iy+1, jx-1), lat_grid(iy+1, jx), lat_grid(iy+1, jx+1))
            min_x = MIN(lon_grid(iy-1, jx-1), lon_grid(iy-1, jx), lon_grid(iy-1, jx+1), &
                        lon_grid(iy, jx-1), lon_grid(iy, jx), lon_grid(iy, jx+1), &
                        lon_grid(iy+1, jx-1), lon_grid(iy+1, jx), lon_grid(iy+1, jx+1))
            max_x = MAX(lon_grid(iy-1, jx-1), lon_grid(iy-1, jx), lon_grid(iy-1, jx+1), &
                        lon_grid(iy, jx-1), lon_grid(iy, jx), lon_grid(iy, jx+1), &
                        lon_grid(iy+1, jx-1), lon_grid(iy+1, jx), lon_grid(iy+1, jx+1))

            if(min_y > max_y) then
                print*, "Error! min y value > max y" 
                stop
            endif
            if(min_x > max_x) then 
                print*, "Error! min x value > max x"
                stop
            endif

            num_loc_counter = 0

            if(max_x - min_x > max_dx) then
                ! print*, "grid cell size >", max_dx, " degrees, likely around lon=0"
                ! print*, "max_x= ", max_x, " min_x = ", min_x
                ! grid_dat(iy/2, jx/2, 1) =  0    
                ! cycle
                Do j = 1, ntot_ims
                    if((lat_grid_ims(j) .ge. min_y) .AND. (lat_grid_ims(j) .le. max_y) .AND. &
                        ((lon_grid_ims(j) .ge. max_x .AND. lon_grid_ims(j) .le. 360.) .OR. &
                        (lon_grid_ims(j) .le. min_x .AND. lon_grid_ims(j) .ge. 0.))) then  
                        num_loc_counter = num_loc_counter + 1
                        grid_dat(iy/2, jx/2, num_loc_counter+1) =  j             
                    endif
                    if (num_loc_counter > num_sub-2) then 
                        print*, "warning more than ",num_sub-2," grids found!"
                        print*, "ymin, ymax, xmin, xmax ", min_y, max_y, min_x, max_x 
                        exit
                    end if
                End do                 
            else               
                Do j = 1, ntot_ims
                    if((lat_grid_ims(j) .ge. min_y) .AND. (lat_grid_ims(j) .le. max_y) .AND. &
                        (lon_grid_ims(j) .ge. min_x) .AND. (lon_grid_ims(j) .le. max_x)) then  
                        num_loc_counter = num_loc_counter + 1
                        grid_dat(iy/2, jx/2, num_loc_counter+1) =  j             
                    endif
                    if (num_loc_counter > num_sub-2) then 
                        print*, "warning more than ",num_sub-2," grids found!"
                        print*, "ymin, ymax, xmin, xmax ", min_y, max_y, min_x, max_x 
                        exit
                    end if
                End do   
            endif         
            grid_dat(iy/2, jx/2, 1) =  num_loc_counter ! first location, num obs
        End do

    End do

    return !grid_dat

End subroutine select_imscells_within_tile_grid
