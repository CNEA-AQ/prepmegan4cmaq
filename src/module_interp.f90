module interpolate_mod

 use netcdf 

 use utils_mod         !utils
 use PROJ_mod          !subroutines for coordinate transformations
 use nc_handler_mod    !functions to deal with netcdf
 implicit none

contains

!!"HOME-MADE" interpolation function:
function interpolate(p,g,inp_file,varname)       result(img2)
 implicit none
 type(grid_type), intent(in)  :: g  !desired grid
 type(proj_type), intent(in)  :: p  !proj of desired grid
 character(*), intent(in)     :: inp_file,varname
 
 real,allocatable :: img2(:,:)      !output array

 integer :: i,j,k

 type(grid_type)  :: GG,GC          !global grid (input grid) &  global grid (CROPPED)
 real,allocatable :: img1(:,:)      !cropped img to be interpolated
 real,allocatable :: lat(:),lon(:)

 integer :: ncid,latid,lonid,varid  !int for netcdf handling

 integer :: is,ie,js,je     !indices that defines subarray
 real    :: px,py,x,y       !dummy variables for temporal coordinates
 real    :: w11,w12,w21,w22 !weights for bilinear interpolation
 real    :: p11,p12,p21,p22 !params. for bilinear interpolation
 real    :: x1,x2,y1,y2     !params. for interpolation
 integer :: i1,i2,j1,j2     !dummy indexes for interpolation

 integer :: scale_x,scale_y

 print*,"  Interpolando: "//trim(inp_file)//"..."

 !Leo inp_file:
 call check(nf90_open(trim(inp_file), nf90_write, ncid ))
   call check( nf90_inq_dimid(ncid, "lat",latid )             )
   call check( nf90_inquire_dimension(ncid, latid, len=GG%ny ))
   call check( nf90_inq_dimid(ncid, "lon",lonid )             )
   call check( nf90_inquire_dimension(ncid, lonid, len=GG%nx ))
   allocate(lat(GG%ny))
   allocate(lon(GG%nx))
   !lat-----------------------------------------------------     
   call check( nf90_inq_varid(ncid,trim("lat"), varid    ))
   !call check( nf90_get_var(ncid, varid , lat(GG%ny:1:-1)))     !lat viene alreves
   call check( nf90_get_var(ncid, varid , lat            ))    
   !lon-----------------------------------------------------     
   call check( nf90_inq_varid(ncid,trim("lon"), varid   ))
   call check( nf90_get_var(ncid, varid , lon           ))
   !------------------------------------------------------
   !Levanto parametros de grilla a interpolar:
   GG%latmin=lat(    1); GG%lonmin=lon(  1)     !lower-left corner?
   GG%latmax=lat(GG%ny); GG%lonmax=lon(GG%nx)   !upper-right corner?

   GG%dy=ABS(GG%latmin-lat(2))  !delta lat
   GG%dx=ABS(GG%lonmin-lon(2))  !delta lon
        !Checkear que sea una grilla regular
        if( ABS(GG%dx - ABS(GG%lonmax-GG%lonmin)/(GG%nx-1)) < 1E-5  ) then; continue;else; print*,"Lon NO es regular.",GG%dx;stop;endif
        if( ABS(GG%dy - ABS(GG%latmax-GG%latmin)/(GG%ny-1)) < 1E-5  ) then; continue;else; print*,"Lat NO es regular.",GG%dy;stop;endif
   !indices de sub-array:
   is=MAX(  1   ,   FLOOR( (g%lonmin-GG%lonmin)/GG%dx) ) !calc min y max indices    
   ie=MIN(GG%nx , CEILING( (g%lonmax-GG%lonmin)/GG%dx) ) !calc min y max indices 
   js=MAX(  1   ,   FLOOR( (g%latmin-GG%latmin)/GG%dy) ) !calc min y max indices
   je=MIN(GG%ny , CEILING( (g%latmax-GG%latmin)/GG%dy) ) !calc min y max indices
   !parametros de grilla:
   GC%nx=ABS(ie-is)+1;  GC%ny=ABS(je-js)+1 
   GC%dx=GG%dx       ;  GC%dy=GG%dy       
   GC%lonmin=lon(is) ;  GC%latmin=lat(js)
   GC%lonmax=lon(ie) ;  GC%latmax=lat(je)

   deallocate(lat)
   deallocate(lon)
   allocate(img1(GC%nx, GC%ny)) 
   !levanto variable a interpolat (es lo que mas tarda)-----     
   call check( nf90_inq_varid(ncid,trim(varname), varid ))
   !call check( nf90_get_var(ncid, varid , img1(1:GC%nx,GC%ny:1:-1), start=[is,GG%ny-je],count=[GC%nx,GC%ny] ) ) !Acá tarda MUCHO..!lat viene alreves
   call check( nf90_get_var(ncid, varid , img1, start=[is,js],count=[GC%nx,GC%ny] ) )
   !--------------------------------------------------------     
 call check(nf90_close(ncid))
 print*,"       File, minval, maxval",trim(inp_file),minval(img1),maxval(img1)
 
 !INTERPOLAR:
 ! Asumo que estoy trabajando con grillas regulares (dx/dy =cte.).    
 ! Asumo que lat y lon estan ordenados de forma creciente.
 
 !Veo si la grilla destino es mas densa o no que la original.
 call xy2ll(p,g%xmin,g%ymin,x1,y1)  !
 call xy2ll(p,g%xmax,g%ymax,x2,y2)  !

 scale_x=CEILING((x2-x1)/(g%nx)/GC%dx)
 scale_y=CEILING((y2-y1)/(g%ny)/GC%dy)

 if (scale_x > 10 .or. scale_y > 10 ) then
    !! Average:
    allocate(img2(g%nx,g%ny))  !array a interpolar:
    !print*,"Average regriding:",scale_x,scale_y
    do i=1,g%nx
       do j=1,g%ny
           px=g%xmin+g%dx*i  !projected coordinate-x
           py=g%ymin+g%dy*j  !projected coordinate-y
                                                                                       
           call xy2ll(p,px,py,x,y)  !I want coordinates in same proj than global file.

           i1=MAX(     1, FLOOR( (x-GC%lonmin) / GC%dx - scale_x*0.5 ) )
           j1=MAX(     1, FLOOR( (y-GC%latmin) / GC%dy - scale_y*0.5 ) )
           i2=MIN( GC%nx, i1+scale_x                                   )
           j2=MIN( GC%ny, j1+scale_y                                   )

           if ( i1 > 1 .and. i2 < GC%nx .and. j1 > 1 .and. j2 < GC%ny ) then
               img2(i,j)=SUM(img1(i1:i2,j1:j2))/((i2-i1)*(j2-j1))  !average
           else
               img2(i,j)=0
           endif
       enddo
    enddo
 else
    !! Bilineal Interp:
    allocate(img2(g%nx,g%ny))  !array a interpolar:
    !print*,"Bilinear interpolation. ",scale_x,scale_y
    do i=1,g%nx
       do j=1,g%ny
           !Position where to interpolate
           px=g%xmin+g%dx*i  !projected coordinate-x
           py=g%ymin+g%dy*j  !projected coordinate-y

           call xy2ll(p,px,py,x,y)  !I want coordinates in same proj than global file.

           !indices:
           i1=FLOOR( (x-GC%lonmin) / GC%dx );  i2=i1+1
           j1=FLOOR( (y-GC%latmin) / GC%dy );  j2=j1+1
           if ( i1 > 1 .and. i2 <= GC%nx .and. j1 > 1 .and. j2 <= GC%ny ) then
              !points (coordinates)    !   p12(i1,j2)    p22(i2,j2)
              x1=GC%lonmin+GC%dx*i1    !       *- - - - - -*            
              x2=GC%lonmin+GC%dx*i2    !       |           |           
              y1=GC%latmin+GC%dy*j1    !       |           |           
              y2=GC%latmin+GC%dy*j2    !       |           |           
              !points (values)         !       |           |           
              p11=img1(i1,j1)          !       *- - - - - -*                                    
              p12=img1(i1,j2)          !   p11(i1,j1)    p21(i2,j1)
              p21=img1(i2,j1)
              p22=img1(i2,j2)
              !weights:
              w11 =(x2 - x )*(y2 - y )/(GC%dx*GC%dy)
              w12 =(x  - x1)*(y2 - y )/(GC%dx*GC%dy)
              w21 =(x2 - x )*(y  - y1)/(GC%dx*GC%dy)
              w22 =(x  - x1)*(y  - y1)/(GC%dx*GC%dy)
              !Bilineal formula:
              img2(i,j)= p11*w11 + p12*w12 + p21*w21 + p22*w22 ! DOT_PRODUCT(p,w)
           else
              img2(i,j)=0.0
           endif
       enddo
    enddo
 endif

 end function

!============================================================================================
!!Funcion de intepolación usando GDAL/OGR..
! function interpolate(p,g,inp_file,varname)   result (img)
!    implicit none
!    type(grid_type) ,intent(in) :: g
!    type(proj_type) ,intent(in) :: p
!    character(*),intent(in) :: inp_file,varname
!    real,allocatable :: img(:,:)
!
!    character(500) :: command
!    character(200) :: proj4
!
!         if ( p%typ == 2 ) then  !Lambert Conformal Conic:
!       proj4="+proj=lcc +lat_1="//trim(rtoa(p%alp))//" +lat_2="//trim(rtoa(p%bet))//" +lon_0="//trim(rtoa(p%gam))//" +lat_0="//trim(rtoa(p%ycent))//" +a=6370000.0 +b=6370000.0 +units=m"
!    else if ( p%typ == 6 ) then  !Polar Secant Stereographic
!       proj4="+proj=stere +lat_ts="//trim(rtoa(p%alp))//" +lon_0="//trim(rtoa(p%gam))//" +a=6370000.0 +b=6370000.0"! +k_0=1.0"
!    else if ( p%typ == 7 ) then  !Equatorial Mercator
!       proj4="+proj=merc +lat_ts="//trim(rtoa(p%alp))//" +lon_0="//trim(rtoa(p%gam))//" +a=6370000.0 +b=6370000.0"
!    else
!       print*, "codigo de proyección invalido.", p%typ; stop
!    end if
!
!    command="gdalwarp -q -overwrite -s_srs 'epsg:4326' -t_srs '"//trim(proj4)//"' -te "//trim(rtoa(g%xmin))//" "//trim(rtoa(g%ymin))//" "//" "//trim(rtoa(g%xmax))//" "//trim(rtoa(g%ymax))//" -tr "//trim(rtoa(g%dx))//" "//trim(rtoa(g%dy))//" -r bilinear -f 'NetCDF' "//trim(inp_file)//" ./tmp.nc "
!
!    print*,"  Interpolando: "//trim(inp_file)//"..."
!    !print*,"$> ",command
!    call system(trim(command))
!
!    allocate(img(g%nx,g%ny))
!    img(:,:)=get2DvarFromNetCDF("./tmp.nc", "Band1", g%nx, g%ny)
! end function
!
!!--------------------------------------------------------------------------------------------

end module

