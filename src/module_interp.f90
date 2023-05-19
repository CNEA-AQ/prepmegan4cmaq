module INTERP_mod
 implicit none

 use netcdf 

 use proj_module
 use netcdf
 use utils_mod         !utils
 use PROJ_mod          !subroutines for coordinate transformations
 use nc_handler_mod    !functions to deal with netcdf
 !use INTERP_mod        !subroutines for interpolation/regridding


 type(grid_type) ,intent(in) :: g
 type(proj_type) ,intent(in) :: p
 character(*),intent(in) :: inpfile,varname
 real,allocatable :: img(:,:)

 integer :: ncid, latid,lonid
 integer :: nr1, nc1, nr2, nc2

 character(*) :: inpfile
 character(16) ::varname

 inpfile='./input/GF3aCrop.nc'
 varname='m20crop'

 call check(nf90_open(trim(inpfile), nf90_write, ncid ))
      
   call check(   nc_inq_dimid(ncid, "lat",latid )  )
   call check(   nc_inq_dimlen(ncid, latid, nr1 )  )
   call check(   nc_inq_dimid(ncid, "lon",lonid )  )
   call check(   nc_inq_dimlen(ncid, lonid, nc1 )  )

   allocate(img(nc1,nr1))
   allocate(lat(nc1,nr1))
   allocate(lon(nc1,nr1))
   !lat----------------------------------------------------------
   call check(   nf90_inq_varid(ncid,trim("lat"), var_id   ))
   call check(   nf90_get_var(ncid, var_id , lat   ))
   !lon----------------------------------------------------------
   call check(   nf90_inq_varid(ncid,trim("lon"), var_id   ))
   call check(   nf90_get_var(ncid, var_id , lon   ))
   !var----------------------------------------------------------
   call check(   nf90_inq_varid(ncid,trim(varname), var_id   ))
   call check(   nf90_get_var(ncid, var_id , img   ))
   !-------------------------------------------------------------

 call check(nf90_close(ncid))

 !(1) Crop file array with boundary box:
 minlat=lat(0,0)        !lower-left corner?
 minlon=lon(0,0)        !lower-left corner?
 maxlat=lat(nr1,nc1)    !upper-right corner?
 maxlon=lon(nr1,nc1)    !upper-right corner?
 dlat=ABS(minlat-lat(1,0))  !ABS(maxlat-minlat)/nr1
 dlon=ABS(minlon-lon(0,1))  !ABS(maxlon-minlon)/nc1



 g%minlat 
 g%minlon
 g%maxlat 
 g%maxlon

 !(2) Create xx1 and yy1 arrays:
 !allocate(xx1(nc1))
 !allocate(yy1(nr1))
 !
 !do i=1,nc1
 !   do j=1, nr1
 !       call ll2xy(p,lon,lat,x,y,)
 !   enddo
 !enddo

 !(3) Interpolate:





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
!
!end module