program test

 use netcdf 

 use proj_mod
 use utils_mod         !utils
 use PROJ_mod          !subroutines for coordinate transformations
 use nc_handler_mod    !functions to deal with netcdf
 use readGRIDDESC_mod

 implicit none

 type(grid_type)  :: grid
 type(proj_type)  :: proj
 character(200)   :: inpfile,varname
 real,allocatable :: img(:,:)
 real,allocatable :: lat(:),lon(:)
 integer :: ny1, nx1

 integer :: ncid, latid,lonid,var_id

 real           :: minlat,minlon,maxlat,maxlon,dlat,dlon
 character(200) :: gridname,griddesc_file

 integer :: is,ie,js,je
 real,allocatable :: img_crop(:,:)
 real,allocatable :: lat_crop(:),lon_crop(:)
 
 !------------------------------------------------------
 !Leo GRIDDESC:
 griddesc_file="GRIDDESC"
 gridname="MERC_TEST"
 call read_GRIDDESC(griddesc_file,gridname, proj, grid)

 call set_additional_proj_params(proj)
 call set_additional_grid_params(proj, grid)

 !------------------------------------------------------
 !Levanto archivo a interpolar:
 inpfile='./input/GF3aCrop.nc'
 varname='m20crop'
 call check(nf90_open(trim(inpfile), nf90_write, ncid ))
   call check( nf90_inq_dimid(ncid, "lat",latid )             )
   call check( nf90_inquire_dimension(ncid, latid, len=ny1 )  )
   call check( nf90_inq_dimid(ncid, "lon",lonid )             )
   call check( nf90_inquire_dimension(ncid, lonid, len=nx1 )  )
   print*,"Global grid: ncols, nrows: ",nx1,ny1
   allocate(img(nx1,ny1))
   allocate(lat(ny1))
   allocate(lon(nx1))
   !lat----------------------------------------------------------
   call check( nf90_inq_varid(ncid,trim("lat"), var_id)   )
   call check( nf90_get_var(ncid, var_id , lat   )        )
   !lon----------------------------------------------------------
   call check( nf90_inq_varid(ncid,trim("lon"), var_id   ))
   call check( nf90_get_var(ncid, var_id , lon   )        )
   !var----------------------------------------------------------
   call check( nf90_inq_varid(ncid,trim(varname), var_id ))
   call check( nf90_get_var(ncid, var_id , img   )        )
   !-------------------------------------------------------------
 call check(nf90_close(ncid))
 !------------------------------------------------------

 !Levanto parametros de grilla a interpolar:
 minlat=lat(ny1)  ; minlon=lon(  1)   !lower-left corner?
 maxlat=lat(  1)  ; maxlon=lon(nx1)   !upper-right corner?
 dlat=ABS(maxlat-lat(2))            !delta lat
 dlon=ABS(minlon-lon(2))            !delta lon

 !Checkear que sea una grilla regular
 if( ABS(dlat - ABS(maxlat-minlat)/ny1) < 0.001  ) then; print*,"Lat OK";else; print*,"Lat NO es regular.";stop;endif
 if( ABS(dlon - ABS(maxlon-minlon)/nx1) < 0.001  ) then; print*,"Lon OK";else; print*,"Lon NO es regular.";stop;endif
 print*,"lat: min max dl",minlat,maxlat,dlat; print*,"lon: min max dl",minlon,maxlon,dlon

 !!(1) Crop file array with boundary box:

 is=MAX(1   ,FLOOR(  ( (grid%lonmin-minlon) )/dlon) )   
 ie=MIN(nx1 ,CEILING(( (grid%lonmax-minlon) )/dlon) )
 je=ny1-MAX(  1 ,FLOOR(  ( (grid%latmin-minlat) )/dlat) ) 
 js=ny1-MIN(ny1 ,CEILING(( (grid%latmax-minlat) )/dlat) )

 print*,"indices for lon",is,ie; print*,"indices for lat",js,je

 grid%nx=ie-is+1
 grid%ny=je-js+1

 allocate(lon_crop(ie-is))
 allocate(lat_crop(je-js))
 allocate(img_crop(ie-is, je-js))

 lon_crop=lon(is:ie)
 lat_crop=lat(js:je)
 img_crop=img(is:ie,js:je)

 print*,"Size of img",shape(img_crop)
 print*,"Shape of grid",grid%nx,grid%ny
 deallocate(img)
 deallocate(lat)
 deallocate(lon)

 ! Create the NetCDF file
 call createNetCDF("tmp.nc",proj,grid,(/'img             '/),(/'units           '/),(/'vardesc                              '/))

 !Abro NetCDF outFile
 call check(nf90_open("tmp.nc", nf90_write, ncid       ))
  call check(nf90_inq_varid(ncid,'img             ',var_id))
  call check(nf90_put_var(ncid, var_id, img_crop ))
 call check(nf90_close( ncid ))
 !Cierro NetCDF outFile



 !(2) Create xx1 and yy1 arrays:
 !allocate(x1(nx1))
 !allocate(y1(ny1))
 !
 !do i=1,nx1
 !   do j=1, ny1
 !       call ll2xy(p,lon,lat,x,y,)
 !   enddo
 !enddo

 !(3) Interpolate:

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
!!        real, intent(in),    allocatable :: Im1(:,:)       !original image
!!        real, intent(in),    allocatable :: X1(:) ,Y1(:)  
!!        real, intent(inout), allocatable :: Im2(:,:)       !interpolated image
!!        real, intent(inout), allocatable :: X2(:) ,Y2(:)
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
!
!end module
