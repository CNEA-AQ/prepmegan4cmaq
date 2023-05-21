module interpolate_mod

 use netcdf 

 use utils_mod         !utils
 use PROJ_mod          !subroutines for coordinate transformations
 use nc_handler_mod    !functions to deal with netcdf
 implicit none

contains

 !Funcion de intepolación usando GDAL/OGR..
 function interpolate(p,g,inp_file,varname)   result (img)
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(*),intent(in) :: inp_file,varname
    real,allocatable :: img(:,:)

    character(500) :: command
    character(200) :: proj4

         if ( p%typ == 2 ) then  !Lambert Conformal Conic:
       proj4="+proj=lcc +lat_1="//trim(rtoa(p%alp))//" +lat_2="//trim(rtoa(p%bet))//" +lon_0="//trim(rtoa(p%gam))//" +lat_0="//trim(rtoa(p%ycent))//" +a=6370000.0 +b=6370000.0 +units=m"
    else if ( p%typ == 6 ) then  !Polar Secant Stereographic
       proj4="+proj=stere +lat_ts="//trim(rtoa(p%alp))//" +lon_0="//trim(rtoa(p%gam))//" +a=6370000.0 +b=6370000.0"! +k_0=1.0"
    else if ( p%typ == 7 ) then  !Equatorial Mercator
       proj4="+proj=merc +lat_ts="//trim(rtoa(p%alp))//" +lon_0="//trim(rtoa(p%gam))//" +a=6370000.0 +b=6370000.0"
    else
       print*, "codigo de proyección invalido.", p%typ; stop
    end if

    command="gdalwarp -q -overwrite -s_srs 'epsg:4326' -t_srs '"//trim(proj4)//"' -te "//trim(rtoa(g%xmin))//" "//trim(rtoa(g%ymin))//" "//" "//trim(rtoa(g%xmax))//" "//trim(rtoa(g%ymax))//" -tr "//trim(rtoa(g%dx))//" "//trim(rtoa(g%dy))//" -r bilinear -f 'NetCDF' "//trim(inp_file)//" ./tmp.nc "

    print*,"  Interpolando: "//trim(inp_file)//"..."
    !print*,"$> ",command
    call system(trim(command))

    allocate(img(g%nx,g%ny))
    img(:,:)=get2DvarFromNetCDF("./tmp.nc", "Band1", g%nx, g%ny)
 end function


!HOME-MADE Interpolation function:
!function interpolate(p,g,inpfile,varname)       result(img2)
! implicit none
! type(grid_type), intent(in)  :: g           !desired grid
! type(proj_type), intent(in)  :: p           !proj of desired grid
! character(*), intent(in)   :: inpfile,varname
! 
! real,allocatable :: img2(:,:)                  !output array
!
! integer :: i,j,k
!
! type(grid_type)  :: GG,GC              !global grid (input grid) &  global grid (CROPPED)
! real,allocatable :: img1(:,:)          !cropped img to be interpolated
! real,allocatable :: lat(:),lon(:)
!
! integer :: ncid,latid,lonid,varid         !int for netcdf handling
!
! integer :: is,ie,js,je                     !usefull indices
! real    :: px,py,x,y
! real    :: w11,w12,w21,w22 !weights for bilinear interpolation
! real    :: p11,p12,p21,p22 !weights for bilinear interpolation
! real    :: x1,x2,y1,y2     !weights for bilinear interpolation
! integer :: i1,i2,j1,j2     !weights for bilinear interpolation
!
! !Leo inpfile:
! call check(nf90_open(trim(inpfile), nf90_write, ncid ))
!   call check( nf90_inq_dimid(ncid, "lat",latid )             )
!   call check( nf90_inquire_dimension(ncid, latid, len=GG%ny ))
!   call check( nf90_inq_dimid(ncid, "lon",lonid )             )
!   call check( nf90_inquire_dimension(ncid, lonid, len=GG%nx ))
!   allocate(lat(GG%ny))
!   allocate(lon(GG%nx))
!   !lat-----------------------------------------------------     
!   call check( nf90_inq_varid(ncid,trim("lat"), varid    ))
!   call check( nf90_get_var(ncid, varid , lat(GG%ny:1:-1)))     !lat viene alreves
!   !lon-----------------------------------------------------     
!   call check( nf90_inq_varid(ncid,trim("lon"), varid   ))
!   call check( nf90_get_var(ncid, varid , lon           ))
!   !------------------------------------------------------
!   !Levanto parametros de grilla a interpolar:
!   GG%latmin=lat(    1); GG%lonmin=lon(  1)     !lower-left corner?
!   GG%latmax=lat(GG%ny); GG%lonmax=lon(GG%nx)   !upper-right corner?
!
!   GG%dy=ABS(GG%latmin-lat(2))            !delta lat
!   GG%dx=ABS(GG%lonmin-lon(2))            !delta lon
!          !Checkear que sea una grilla regular
!          if( ABS(GG%dy - ABS(GG%latmax-GG%latmin)/GG%ny) < 1E-5  ) then; print*,"Lat OK";else; print*,"Lat NO es regular.";stop;endif
!          if( ABS(GG%dx - ABS(GG%lonmax-GG%lonmin)/GG%nx) < 1E-5  ) then; print*,"Lon OK";else; print*,"Lon NO es regular.";stop;endif
!   !------------------------------------------------------
!   !!(1) Genero parametros de grilla recortada:
!   is= MAX(1     ,  FLOOR(( (g%lonmin-GG%lonmin) )/GG%dx) ) !calc min y max indices    
!   ie= MIN(GG%nx ,CEILING(( (g%lonmax-GG%lonmin) )/GG%dx) ) !calc min y max indices 
!   js= MAX(  1   ,  FLOOR(( (g%latmin-GG%latmin) )/GG%dy) ) !calc min y max indices
!   je= MIN(GG%ny ,CEILING(( (g%latmax-GG%latmin) )/GG%dy) ) !calc min y max indices
!   
!   GC%nx=ABS(ie-is)+1;  GC%ny=ABS(je-js)+1 
!   GC%dx=GG%dx       ;  GC%dy=GG%dy       
!   GC%lonmin=lon(is) ;  GC%latmin=lat(js)
!   GC%lonmax=lon(ie) ;  GC%latmax=lat(je)
!
!   deallocate(lat)
!   deallocate(lon)
!       print*,"GC: indices for lon",is,"-",ie;   
!       print*,"GC: indices for lat",js,"-",je
!       !print*,"GG: ncols, nrows:  ",GG%nx,GG%ny
!       !print*,"GG: lat: min max dl",GG%latmin,GG%latmax,GG%dy; !print*,"dl, dl",GG%dx, ABS(GG%lonmax-GG%lonmin)/(GG%nx)
!       !print*,"GG: lon: min max dl",GG%lonmin,GG%lonmax,GG%dx; !print*,"dl, dl",GG%dy, ABS(GG%latmax-GG%latmin)/(GG%ny)
!       print*,"GC: ncols, nrows:  ",GC%nx,GC%ny
!       print*,"GC: lat: min max dl",GC%latmin,GC%latmax,GC%dy 
!       print*,"GC: lon: min max dl",GC%lonmin,GC%lonmax,GC%dx
!       print*,"g : ncols, nrows:  ",g%nx,g%ny 
!       print*,"g : lat: min max dl",g%latmin,g%latmax,g%dy 
!       print*,"g : lon: min max dl",g%lonmin,g%lonmax,g%dx
!   !------------------------------------------------------
!   !!(2) Crop Img original                               
!   !allocate(img1(GC%nx, GC%ny))
!   !img1=img(is:ie,je:js:-1)          !aca invierto el orden de "y" para que quede ordenado de forma creciente.                      
!   allocate(img1(GC%nx, GC%ny)) 
!
!   !var (esto es lo que mas tarda)--------------------------     
!   call check( nf90_inq_varid(ncid,trim(varname), varid ))
!   call check( nf90_get_var(ncid, varid , img1(1:GC%nx,GC%ny:1:-1), start=[is,GG%ny-je],count=[GC%nx,GC%ny] ) ) !Acá tarda MUCHO..
!   !--------------------------------------------------------     
! call check(nf90_close(ncid))
! 
!     !!Esto es solo de testeo de grilla cropeada
!     !! Create the NetCDF file
!     !call createNetCDF("tmp_GC.nc",p,GC,(/'img             '/),(/'units           '/),(/'vardesc                              '/))
!     !!Abro NetCDF outFile
!     !call check(nf90_open("tmp_GC.nc", nf90_write, ncid       ))
!     ! call check(nf90_inq_varid(ncid,'img             ',varid))
!     ! call check(nf90_put_var(ncid, varid, img1 ))
!     !call check(nf90_close( ncid ))
!     !!Cierro NetCDF outFile
!
! !INTERPOLAR:
! ! Asumo que estoy trabajando con grillas regulares (dx/dy =cte.).    
! ! Asumo que lat y lon estan ordenados de forma creciente.
! 
! allocate(img2(g%nx,g%ny))  !array a interpolar:
! do i=1,g%nx
!    do j=1,g%ny
!        
!        !Position where to interpolate
!        px=g%xmin+g%dx*i  !projected coordinate-x
!        py=g%ymin+g%dy*j  !projected coordinate-y
!
!        call xy2ll(p,px,py,x,y)  !I want coordinates in same proj than global file.
!        !print*,"px,py,x,y: ",px,py,x,y
!        
!        !indices:
!        i1=FLOOR( (x-GC%lonmin) / GC%dx );  i2=i1+1
!        j1=FLOOR( (y-GC%latmin) / GC%dy );  j2=j1+1
!        if ( i1 > 1 .and. i2 <= GC%nx .and. j1 > 1 .and. j2 <= GC%ny ) then
!
!            !points (coordinates)    !   p12(i1,j2)    p22(i2,j2)
!            x1=GC%lonmin+GC%dx*i1    !       *- - - - - -*            
!            x2=GC%lonmin+GC%dx*i2    !       |           |           
!            y1=GC%latmin+GC%dy*j1    !       |           |           
!            y2=GC%latmin+GC%dy*j2    !       |           |           
!            !points (values)         !       |           |           
!            p11=img1(i1,j1)          !       *- - - - - -*                                    
!            p12=img1(i1,j2)          !   p11(i1,j1)    p21(i2,j1)
!            p21=img1(i2,j1)
!            p22=img1(i2,j2)
!            !weights:
!            w11 =(x2 - x )*(y2 - y )/(GC%dx*GC%dy)
!            w12 =(x  - x1)*(y2 - y )/(GC%dx*GC%dy)
!            w21 =(x2 - x )*(y  - y1)/(GC%dx*GC%dy)
!            w22 =(x  - x1)*(y  - y1)/(GC%dx*GC%dy)
!                !if (p11 > 0 .or. p21 > 0) then
!                !    print*,"i1,i2,i1,i2:",i1,i2,i1,i2
!                !    print*,"x,y,x1,y1,x2,y2:",x,y,x1,x2,y1,y2
!                !    print*,"p11,p12,p21,p22:",p11,p12,p21,p22
!                !    print*,"w11,w12,w21,w22:",w11,w12,w21,w22
!                !endif
!            !Bilineal formula:
!            img2(i,j)= p11*w11 + p12*w12 + p21*w21 + p22*w22 ! DOT_PRODUCT(p,w)
!        else
!            img2(i,j)=0.0
!        endif
!    enddo
! enddo
!
! end function


end module


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
!
!
!
!end function
!
!
!subroutine interpolate(Im1, X1, Y1, Im2, X2, Y2)
!  implicit none
!  real, intent(in) :: X1(:), Y1(:), X2(:), Y2(:)
!  real, intent(in)    :: Im1(:,:)
!  real, intent(inout) :: Im2(:,:)
!  integer :: N1, N2, i, j
!  real :: u, v
!  real :: w1, w2, w3, w4
!  
!  N1 = size(Im1, 1) ! Number of rows in Im1
!  N2 = size(Im2, 1) ! Number of rows in Im2
!  
!  ! Perform bilinear interpolation
!  do i = 1, N2
!    do j = 1, N2
!      where(X1 >= 1.0_r8 .and. X1 <= N1 .and. Y1 >= 1.0_r8 .and. Y1 <= N1)
!        ! Calculate integer indices and interpolation weights
!        u = X1(j) - floor(X1(j))
!        v = Y1(i) - floor(Y1(i))
!        
!        ! Perform bilinear interpolation
!        w1 = (1.0_r8 - u) * (1.0_r8 - v)
!        w2 = u * (1.0_r8 - v)
!        w3 = (1.0_r8 - u) * v
!        w4 = u * v
!        
!        Im2(i,j) = w1 * Im1(int(X1(j)), int(Y1(i))) &
!                 + w2 * Im1(int(X1(j))+1, int(Y1(i))) &
!                 + w3 * Im1(int(X1(j)), int(Y1(i))+1) &
!                 + w4 * Im1(int(X1(j))+1, int(Y1(i))+1)
!      elsewhere
!        ! Set Im2 to a default value when (X1, Y1) is outside the valid range
!        Im2(i,j) = 0
!      end where
!    end do
!  end do
!end subroutine
!
!
!!        !Bilinear Interpolation:
!!        !
!!        !  Im1(X1,Y1) --> Im2(X2,Y2)
!!        !
!!        !
!!        !
!!        !
!!
!!
!!subroutine interpolate(Im1,X1,Y1,Im2,X2,Y2)
!!        implicit none
!!        real, intent(in),    allocatable :: Im1(:), X1(:) ,Y1(:)   !original image
!!        real, intent(inout), allocatable :: Im2(:), X2(:) ,Y2(:)   !interpolated image
!!        
!!        x1 =   FLOOR(x)
!!        x2 = CEILING(x)
!!        y1 =   FLOOR(y)
!!        y2 = CEILING(y)
!!                                                               
!!        w(1) =  (x2-x )*(y2-y )/(x2-x1)/(y2-y1)
!!        w(2) =  (x -x1)*(y2-y )/(x2-x1)/(y2-y1)
!!        w(3) =  (x2-x )*(y2-y )/(x2-x1)/(y2-y1)
!!        w(4) =  (x2-x )*(y2-y )/(x2-x1)/(y2-y1)
!!
!!        Im2=dot_product(Im1,w)
!!
!!        wx1 = x2-x
!!        wy1 = y2-y
!!                                                                       
!!        wx2 = 1 - wx1
!!        wy2 = 1 - wy1
!!                                                                       
!!        FORALL( jj=1 : NPts)
!!                                                                       
!!           Im2(jj) = wy1(jj)*(
!!                              wx1(jj)*Im1(y1(jj),x1(jj))  + wx2(jj)*
!!           $                          Im1(y1(jj),x2(jj))) + wy2(jj)*
!!                             (
!!                             wx1(jj)*
!!           $                          Im1(y2(jj),x1(jj))  + wx2(jj)*
!!                                      Im1(y2(jj),x2(jj))
!!                                      )
!!                                                                       
!!        END FORALL
!!
!!
!!subroutine
!!
!!SUBROUTINE L2DINTERPOL(IntIm,Image,x,y,NPts,M,N)
!!  implicit none
!!
!!  mwSize, PARAMETER             :: dp = kind(0.d0) ! Double precision
!!  mwSize                        :: NPts, M,N       ! Input
!!  REAL(dp),DIMENSION(Npts)      :: x,y             ! Input
!!  REAL(dp),DIMENSION(M,N)       :: Image           ! Input
!!  mwSize                        :: jj
!!  mwSize, DIMENSION(NPts)       :: x1,y1,x2,y2      
!!  REAL(dp),DIMENSION(Npts)      :: wx1,wx2,wy1,wy2
!!  REAL(dp),DIMENSION(Npts)      :: IntIm           ! Output
!!
!!
!!  x1 = FLOOR(x)
!!  x2 = CEILING(x)
!!  y1 = FLOOR(y)
!!  y2 = CEILING(y)
!!
!!  wx1 = x2-x
!!  wy1 = y2-y
!!
!!  wx2 = 1 - wx1
!!  wy2 = 1 - wy1
!!
!!  FORALL( jj=1 : NPts)
!!
!!     IntIm(jj) = wy1(jj)*(wx1(jj)*Image(y1(jj),x1(jj))+wx2(jj)*
!! $        Image(y1(jj),x2(jj))) + wy2(jj)*(wx1(jj)*
!!     $        Image(y2(jj),x1(jj))+wx2(jj)*Image(y2(jj),x2(jj)))
!!
!!  END FORALL
!!END 
!!
