program test

 use netcdf 
 use utils_mod         !utils
 use PROJ_mod          !subroutines for coordinate transformations
 use nc_handler_mod    !functions to deal with netcdf
 use readGRIDDESC_mod

 implicit none
 integer :: i,j,k
 type(grid_type)  :: GG,GC  !global grid (input grid) &  global grid (CROPPED)
 type(grid_type)  :: grid   !desired grid
 type(proj_type)  :: proj   !proj of desired grid
 character(200)   :: inpfile,varname
 real,allocatable :: img2(:,:)
 real,allocatable :: img1(:,:)          !cropped img to be interpolated
 real,allocatable :: lat(:),lon(:)

 character(200) :: gridname, griddesc_file  !for GRIDDESC reading
 integer :: ncid,latid,lonid,var_id         !int for netcdf handling
 real :: x,y,ln,lt

 integer :: is,ie,js,je                     !usefull indices
 real :: px,py
 real    :: w11,w12,w21,w22 !weights for bilinear interpolation
 real    :: p11,p12,p21,p22 !weights for bilinear interpolation
 real    :: x1,x2,y1,y2     !weights for bilinear interpolation
 integer :: i1,i2,j1,j2     !weights for bilinear interpolation
        
 !------------------------------------------------------
 !Leo GRIDDESC:
 griddesc_file="GRIDDESC"
 gridname="MERC_TEST"
 call read_GRIDDESC(griddesc_file,gridname, proj, grid)
 call set_additional_proj_params(proj)
 call set_additional_grid_params(proj, grid)
 !------------------------------------------------------
 print*,"(0) Leo grilla a grilla a interpolar -----------------"
 !Levanto archivo a interpolar:
 inpfile='./input/GF3aCrop.nc'
 varname='m20crop'
 call check(nf90_open(trim(inpfile), nf90_write, ncid ))
   call check( nf90_inq_dimid(ncid, "lat",latid )             )
   call check( nf90_inquire_dimension(ncid, latid, len=GG%ny ))
   call check( nf90_inq_dimid(ncid, "lon",lonid )             )
   call check( nf90_inquire_dimension(ncid, lonid, len=GG%nx ))
   allocate(lat(GG%ny))
   allocate(lon(GG%nx))
   !lat-----------------------------------------------------     
   call check( nf90_inq_varid(ncid,trim("lat"), var_id)   )
   call check( nf90_get_var(ncid, var_id , lat   )        )
   !lon-----------------------------------------------------     
   call check( nf90_inq_varid(ncid,trim("lon"), var_id   ))
   call check( nf90_get_var(ncid, var_id , lon   )        )
   !------------------------------------------------------
   !Levanto parametros de grilla a interpolar:
   GG%latmin=lat(GG%ny); GG%lonmin=lon(  1)     !lower-left corner?
   GG%latmax=lat(  1)  ; GG%lonmax=lon(GG%nx)   !upper-right corner?
   GG%dy=ABS(GG%latmax-lat(2))            !delta lat
   GG%dx=ABS(GG%lonmin-lon(2))            !delta lon
          !Checkear que sea una grilla regular
          if( ABS(GG%dy - ABS(GG%latmax-GG%latmin)/GG%ny) < 1E-5  ) then; print*,"Lat OK";else; print*,"Lat NO es regular.";stop;endif
          if( ABS(GG%dx - ABS(GG%lonmax-GG%lonmin)/GG%nx) < 1E-5  ) then; print*,"Lon OK";else; print*,"Lon NO es regular.";stop;endif
                !print*,"GG: ncols, nrows:  ",GG%nx,GG%ny
                !print*,"GG: lat: min max dl",GG%latmin,GG%latmax,GG%dy; !print*,"dl, dl",GG%dx, ABS(GG%lonmax-GG%lonmin)/(GG%nx)
                !print*,"GG: lon: min max dl",GG%lonmin,GG%lonmax,GG%dx; !print*,"dl, dl",GG%dy, ABS(GG%latmax-GG%latmin)/(GG%ny)
   !------------------------------------------------------
   !!(1) Genero parametros de grilla recortada:
   is=MAX(1     ,  FLOOR(( (grid%lonmin-GG%lonmin) )/GG%dx) ) !calc min y max indices    
   ie=MIN(GG%nx ,CEILING(( (grid%lonmax-GG%lonmin) )/GG%dx) ) !calc min y max indices 
   js=MAX(  1   ,  FLOOR(( (grid%latmin-GG%latmin) )/GG%dy) ) !calc min y max indices
   je=MIN(GG%ny ,CEILING(( (grid%latmax-GG%latmin) )/GG%dy) ) !calc min y max indices
        print*,"GC: indices for lon",is,"-",ie; 
        print*,"GC: indices for lat",js,"-",je
   GC%nx=ABS(ie-is)+1;  GC%ny=ABS(je-js)+1                                                                                              
   GC%dx=GG%dx       ;  GC%dy=GG%dy                                                                                                     
   GC%lonmin=lon(is)                                                                                                                    
   GC%lonmax=lon(ie)                                                                                                                    
   GC%latmin=lat(GG%ny-js)        !esto solo por que lat está alreves                                                               
   GC%latmax=lat(GG%ny-je)        !esto solo por que lat está alreves                                                               
   deallocate(lat)                                                                                                                        
   deallocate(lon)                                                                                                                      
       print*,"GC: ncols, nrows:  ",GC%nx,GC%ny                                                                                             
       print*,"GC: lat: min max dl",GC%latmin,GC%latmax,GC%dy;                                                                              
       print*,"GC: lon: min max dl",GC%lonmin,GC%lonmax,GC%dx                                                                               
       print*,"g : ncols, nrows:  ",grid%nx,grid%ny                                                                                         
       print*,"g : lat: min max dl",grid%latmin,grid%latmax,grid%dy;                                                                        
       print*,"g : lon: min max dl",grid%lonmin,grid%lonmax,grid%dx                                                                         
   !------------------------------------------------------                                                                              
   !!(2) Crop Img original                                                                                                               
   !print*,"(2) Crop Img original ---------------------------------"                                                                     
   !allocate(img1(GC%nx, GC%ny))                                                                                                     
   !img1=img(is:ie,je:js:-1)          !aca invierto el orden de "y" para que quede ordenado de forma creciente.                      
   allocate(img1(GC%nx, GC%ny))
   
   !var (esto es lo que mas tarda)--------------------------     
   call check( nf90_inq_varid(ncid,trim(varname), var_id ))
   call check( nf90_get_var(ncid, var_id , img1, start=[is,js],count=[GC%nx,GC%ny] ) )        !Acá tarda MUCHO..
   !--------------------------------------------------------     
 call check(nf90_close(ncid))
 
     !Esto es solo de testeo de grilla cropeada
     ! Create the NetCDF file
     call createNetCDF("tmp_GC.nc",proj,GC,(/'img             '/),(/'units           '/),(/'vardesc                              '/))
     !Abro NetCDF outFile
     call check(nf90_open("tmp_GC.nc", nf90_write, ncid       ))
      call check(nf90_inq_varid(ncid,'img             ',var_id))
      call check(nf90_put_var(ncid, var_id, img1 ))
     call check(nf90_close( ncid ))
     !Cierro NetCDF outFile

!INTERPOLAR:
! Asumo que estoy trabajando con grillas regulares (dx/dy =cte.).    
! Asumo que lat y lon estan ordenados de forma creciente.
 
 allocate(img2(grid%nx,grid%ny))  !array a interpolar:
 do i=1,grid%nx
        do j=1,grid%ny
        
        !Position where to interpolate
        px=grid%xmin+grid%dx*i  !projected coordinate-x
        py=grid%ymin+grid%dy*j  !projected coordinate-y

        call xy2ll(proj,px,py,x,y)  !I want coordinates in same proj than global file.
        !print*,"px,py,x,y: ",px,py,x,y
        
        !indices:
        i1=FLOOR( (x-GC%lonmin) / GC%dx );  i2=i1+1
        j1=FLOOR( (y-GC%latmin) / GC%dy );  j2=j1+1
        if ( i1 > 1 .and. i2 <= GC%nx .and. j1 > 1 .and. j2 <= GC%ny ) then

            !points (coordinates)        !   p12(i1,j2)    p22(i2,j2)
            x1=GC%lonmin+GC%dx*i1        !       *- - - - - -*            
            x2=GC%lonmin+GC%dx*i2        !       |           |           
            y1=GC%latmin+GC%dy*j1        !       |           |           
            y2=GC%latmin+GC%dy*j2        !       |           |           
            !points (values)             !       |           |           
            p11=img1(i1,j1)          !       *- - - - - -*                                    
            p12=img1(i1,j2)          !   p11(i1,j1)    p21(i2,j1)
            p21=img1(i2,j1)
            p22=img1(i2,j2)
            !weights:
            w11 =(x2 - x )*(y2 - y )/(GC%dx*GC%dy)  !(x2 - x1)/(y2 - y1) !
            w12 =(x  - x1)*(y2 - y )/(GC%dx*GC%dy)  !(x2 - x1)/(y2 - y1) !
            w21 =(x2 - x )*(y  - y1)/(GC%dx*GC%dy)  !(x2 - x1)/(y2 - y1) !
            w22 =(x  - x1)*(y  - y1)/(GC%dx*GC%dy)  !(x2 - x1)/(y2 - y1) !
                if (p11 > 0 .or. p21 > 0) then
                    print*,"i1,i2,i1,i2:",i1,i2,i1,i2
                    print*,"x,y,x1,y1,x2,y2:",x,y,x1,x2,y1,y2
                    print*,"p11,p12,p21,p22:",p11,p12,p21,p22
                    print*,"w11,w12,w21,w22:",w11,w12,w21,w22
                endif

            !Bilineal formula:
            img2(i,j)= p11*w11 + p12*w12 + p21*w21 + p22*w22 ! DOT_PRODUCT(p,w)
        else
            img2(i,j)=0.0
        endif
        !where(X1 >= 1.0_r8 .and. X1 <= N1 .and. Y1 >= 1.0_r8 .and. Y1 <= N1)
        !  ! Calculate integer indices and interpolation weights
        !  u = X1(j) - floor(X1(j))
        !  v = Y1(i) - floor(Y1(i))
        !  
        !  ! Perform bilinear interpolation
        !  w1 = (1.0_r8 - u) * (1.0_r8 - v)
        !  w2 = u * (1.0_r8 - v)
        !  w3 = (1.0_r8 - u) * v
        !  w4 = u * v
        !  
        !  Im2(i,j) = w1 * Im1(int(X1(j)), int(Y1(i))) &
        !           + w2 * Im1(int(X1(j))+1, int(Y1(i))) &
        !           + w3 * Im1(int(X1(j)), int(Y1(i))+1) &
        !           + w4 * Im1(int(X1(j))+1, int(Y1(i))+1)
        !elsewhere
        !  ! Set Im2 to a default value when (X1, Y1) is outside the valid range
        !  Im2(i,j) = 0
        !end where
        enddo
 enddo
!Esto es solo de testeo:
! Create the NetCDF file
call createNetCDF("tmp.nc",proj,grid,(/'img             '/),(/'units           '/),(/'vardesc                              '/))
!Abro NetCDF outFile
call check(nf90_open("tmp.nc", nf90_write, ncid       ))
 call check(nf90_inq_varid(ncid,'img             ',var_id))
 call check(nf90_put_var(ncid, var_id, img2 ))
call check(nf90_close( ncid ))
!Cierro NetCDF outFile



end program


!contains
!
!function interpolate(p,g,inpfile)   result (img)
!  implicit none
!  type(grid_type) ,intent(in) :: g
!  type(proj_type) ,intent(in) :: p
!  character(*),intent(in) :: inpfile,varname
!  real,allocatable :: img(:,:)
!
!
!
!end function
!end module
