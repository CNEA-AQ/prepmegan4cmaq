!------------------------------------------------------
! Based on UCI-BAI-MEGAN code and megan_bio_emis.f90
! 
! author: Ramiro A. Espada. April 2023.
!------------------------------------------------------
program prepmegan4cmaq

  use netcdf
  implicit none

  type proj_type
     character(16)    :: pName     ! nombre de la proyeccion
     character(7)     :: typ_str   ! String  code for projection TYPE
     integer          :: typ       ! Integer code for projection TYPE
     real             :: ref_lat,ref_lon,truelat1,truelat2,stand_lon,pole_lat,pole_lon
     character(125)   :: proj4     ! PROJ4 srs definition.
  end type proj_type

  type grid_type
      character(12)    :: gName        !grid-name
      integer          :: nx,ny,nz     !number of cells in x-y direction (ncols, nrows, nlevs)
      real             :: dx,dy        !x-y cell dimension (x_cell, y_cell)
      real             :: xmin,ymin,xmax,ymax,xc,yc
      real             :: lonmin,latmin,lonmax,latmax
   end type grid_type

  type(proj_type) :: proj
  type(grid_type) :: grid

  integer :: status,iostat
  integer :: i,j,k
  integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,var_id
  logical :: file_exists

  character(len=17)  ::start_date,end_date
  character(200) :: griddesc_file, crop_frac_file,grass_frac_file,shrub_frac_file,tree_frac_file,nl_tree_frac_file,bl_tree_frac_file,tp_tree_frac_file,ecotype_file,laiv_files

  namelist/control/start_date,end_date,griddesc_file, crop_frac_file,grass_frac_file,shrub_frac_file,tree_frac_file,nl_tree_frac_file, bl_tree_frac_file,tp_tree_frac_file,ecotype_file,laiv_files!,lai_file,veg_cov_file,wrfout_file

  !Leo namelist:
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'prepmegan4cmaq: failed to read namelist; error = ',iostat
    stop
  end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file, proj, grid)    !(!) TO-DO: mejorar esta funcion basado en lo que haga IOAPI

  !`MEGAN_CTS` (*Canopy Type Fractions*) 
  call build_CT3(grid,proj,crop_frac_file,tree_frac_file,grass_frac_file,shrub_frac_file,nl_tree_frac_file,tp_tree_frac_file)
  
  !`MEGAN_LAI` (Leaf Area Index).
  !call build_LAI(grid,proj,lai_file)
  
  !`MEGAN_LDF` (*Light Dependence Fractions*)
 

  !`MEGAN_EFS` (emission factors).


print*, "========================================="
print*, " prepmegan4cmaq: Completed successfully"
print*, "========================================="

contains

 subroutine check(status)
   integer, intent(in) :: status
   if (status /= nf90_noerr) then
     write(*,*) nf90_strerror(status)
     stop 'netcdf error'
   end if
 end subroutine check

 !Interfaz a "date"
 function date(date_str, fmt_str) result(output)
   implicit none
   character(*), intent(in) :: date_str, fmt_str
   character(256)           :: command
   character(20)            :: output
   command="date -d "//trim(date_str)//" '+"//trim(fmt_str)//"'  > tmp_date.txt"
   call system( trim(command) )
   !print*,trim(command)
   open(9, file='tmp_date.txt', status='old',action='read'); read(9, '(A)', iostat=status) output;  close(9)
   call system('rm tmp_date.txt')
 end function

 function atoi(str)     !string -> int
   implicit none
   character(len=*), intent(in) :: str
   integer :: atoi
   read(str,*) atoi
 end function
 function itoa(i)       !int -> string
    implicit none
    integer, intent(in) :: i
    character(len=20) :: itoa
    write(itoa, '(i0)') i
    itoa = adjustl(itoa)
 end function
 function rtoa(r)       !real -> string
    implicit none
    real, intent(in) :: r
    character(len=16) :: rtoa
    write(rtoa, '(F16.3)') r
    rtoa = adjustl(rtoa)
 end function

 subroutine read_GRIDDESC(griddescFile, p, g)
        implicit none
        character(200),intent(in) :: griddescFile
        type(proj_type) ,intent(inout) :: p
        type(grid_type) ,intent(inout) :: g
        character(10) :: dummyvar
        open(2,file=griddescFile, status='old', action='read')                                  !GRIDDESC:
           read(2,*) dummyvar;                                                   !' '
           read(2,*) p%pName;                                                    !projName
           read(2,*) p%typ,p%truelat1,p%truelat2,p%stand_lon,p%ref_lon,p%ref_lat !map_proj truelat1 truelat2 stand_lon ref_lon ref_lat
           read(2,*) dummyvar;                                                   !' '
           read(2,*) g%gName;                                                    !gridName
           read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%nx,g%ny                   !projName xorig yorig xcell ycell nrows ncols
        close(2)
 end subroutine

subroutine build_CT3(g,p,cropfile,treefile,grassfile,shrubfile,nltreefile,troptreefile)
   implicit none
   type(grid_type) ,intent(in) :: g
   type(proj_type) ,intent(in) :: p
   character(len=200) :: cropfile,treefile,grassfile,shrubfile,nltreefile,troptreefile
   !real, allocatable :: crop(:,:),shrub(:,:),grass(:,:),tree(:,:),nl_tree(:,:),trop_tree(:,:),bl_tree(:,:),laiv(:,:)
   real, allocatable :: CTS(:,:,:)
   integer :: date_time_dim_id,var_dim_id,lev_id,tstep_dim_id,var_id,ncid
   character(len=6) :: outfile
   character(len=25),allocatable :: var_list(:),var_desc(:)
   character(len=200),allocatable :: file_list(:)
   integer :: nvars

   outfile='CT3.nc'
   nvars=6

   allocate(var_list(nvars))  
   allocate(var_desc(nvars))  
   allocate(file_list(nvars))  
   allocate(CTS(grid%nx,grid%ny,nvars))  

   file_list=(/treefile,shrubfile,cropfile, grassfile, nltreefile, troptreefile/)      !No cambiar, EL ORDEN IMPORTA (do not touch, the order matters)!
   
   !Levanto netcdf input files
   do k=1,nvars
        print*,"k=",k
       call check(nf90_open( trim(file_list(k)), nf90_write, ncid   ))
       call check(   nf90_inq_varid(ncid,'Band1', var_id     ))
       call check(   nf90_get_var(ncid,  var_id , CTS(:,:,k)  ))
       call check(nf90_close(ncid                             ))
   enddo

   var_list =(/ "shrub_frac             ","crop_frac              ","grass_frac             ","nl_tree_frac           ","trop_tree_frac         ","bl_tree_frac           "/)
   var_desc =(/ "shrub fraction         ","crop fraction          ","grass fraction         ","needle tree fraction   ","tropical tree fraction ","broadleaf tree fraction"/) 
   ! Create the NetCDF file
   call check(nf90_create(outFile, NF90_CLOBBER, ncid))
       ! Defino dimensiones
       call check(nf90_def_dim(ncid, "TSTEP"    ,   1    , tstep_dim_id    ))
       call check(nf90_def_dim(ncid, "DATE_TIME",   2    , date_time_dim_id))
       call check(nf90_def_dim(ncid, "COL"      , grid%nx, col_dim_id      ))
       call check(nf90_def_dim(ncid, "ROW"      , grid%ny, row_dim_id      ))
       call check(nf90_def_dim(ncid, "LAY"      ,   1    , lay_dim_id      ))
       call check(nf90_def_dim(ncid, "VAR"      , nvars  , var_dim_id      ))

       !Defino variables
       call check(nf90_def_var(ncid,"TFLAG",NF90_FLOAT    , [date_time_dim_id,var_dim_id,tstep_dim_id], var_id))
       call check(nf90_put_att(ncid, var_id, "units"      , "<YYYYDDD,HHMMSS>" ))
       call check(nf90_put_att(ncid, var_id, "long_name"  , "TFLAG           " ))
       call check(nf90_put_att(ncid, var_id, "var_desc"   , "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                "))

       do k=1, nvars
         call check(nf90_def_var(ncid, trim(var_list(k)) ,NF90_FLOAT , [col_dim_id,row_dim_id,lay_dim_id,tstep_dim_id], var_id)) !
         call check(nf90_put_att(ncid, var_id,"units"    , "nondimension    " ))
         call check(nf90_put_att(ncid, var_id,"long_name", trim(var_list(k))  ))
         call check(nf90_put_att(ncid, var_id,"var_desc" , trim(var_desc(k))  ))
       end do
    call check(nf90_enddef(ncid))
    !End NetCDF define mode

    !Abro NetCDF outFile
    call check(nf90_open(outFile, nf90_write, ncid       ))
       !shrub
       call check(nf90_inq_varid(ncid,"shrub_frac"     ,var_id))
       call check(nf90_put_var(ncid, var_id, CTS(:,:,2)  ))        
       !crop
       call check(nf90_inq_varid(ncid,"crop_frac"     ,var_id))
       call check(nf90_put_var(ncid, var_id, CTS(:,:,3)  ))        
       !grass
       call check(nf90_inq_varid(ncid,"grass_frac"     ,var_id))
       call check(nf90_put_var(ncid, var_id, CTS(:,:,4)  ))        
       !needleleaf tree
       call check(nf90_inq_varid(ncid,"nl_tree_frac"   ,var_id))
       call check(nf90_put_var(ncid, var_id, CTS(:,:,1) * (1-CTS(:,:,6)) * CTS(:,:,5)     )) 
       !boradleaf tree
       call check(nf90_inq_varid(ncid,"bl_tree_frac"   ,var_id))               
       call check(nf90_put_var(ncid, var_id, CTS(:,:,1) * (1-CTS(:,:,6)) * (1-CTS(:,:,5)) ))
       !tropical tree
       call check(nf90_inq_varid(ncid,"trop_tree_frac" ,var_id))
       call check(nf90_put_var(ncid, var_id, CTS(:,:,1)*CTS(:,:,6)                    )) 
       !TFLAG:
       call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
       call check(nf90_put_var(ncid, var_id, (/0000000,000000 /) ))

       !! Defino attributos
       !call check(nf90_put_att(ncid, nf90_global,"IOAPI_VERSION", "ioapi-3.2: \$Id: init3" ))
       !call check(nf90_put_att(ncid, nf90_global,"EXEC_ID"  , "????????????????"   ))
       !call check(nf90_put_att(ncid, nf90_global,"FTYPE"    , 1                    ))
       !call check(nf90_put_att(ncid, nf90_global,"SDATE"    , atoi(YYYY//DDD)      ))   !int
       !call check(nf90_put_att(ncid, nf90_global,"STIME"    , 000000               ))
       !call check(nf90_put_att(ncid, nf90_global,"WDATE"    , 2023001              ))
       !call check(nf90_put_att(ncid, nf90_global,"WTIME"    , 000000               ))
       !call check(nf90_put_att(ncid, nf90_global,"CDATE"    , 2023001              ))
       !call check(nf90_put_att(ncid, nf90_global,"CTIME"    , 000000               ))
       !call check(nf90_put_att(ncid, nf90_global,"TSTEP"    , 10000                ))
       !call check(nf90_put_att(ncid, nf90_global,"NTHIK"    , 1                    ))
       !call check(nf90_put_att(ncid, nf90_global,"NCOLS"    , grid%nx              ))
       !call check(nf90_put_att(ncid, nf90_global,"NROWS"    , grid%ny              ))
       !call check(nf90_put_att(ncid, nf90_global,"NLAYS"    , 1                    ))!grid%nz
       !call check(nf90_put_att(ncid, nf90_global,"NVARS"    , nvars                ))
       !call check(nf90_put_att(ncid, nf90_global,"GDTYP"    , 1                    ))
       !call check(nf90_put_att(ncid, nf90_global,"P_ALP"    , -50.                 ))
       !call check(nf90_put_att(ncid, nf90_global,"P_BET"    , -20.                 ))
       !call check(nf90_put_att(ncid, nf90_global,"P_GAM"    , -65.                 ))
       !call check(nf90_put_att(ncid, nf90_global,"XCENT"    , proj%ref_lon         ))
       !call check(nf90_put_att(ncid, nf90_global,"YCENT"    , proj%ref_lat         ))
       !call check(nf90_put_att(ncid, nf90_global,"XORIG"    , grid%xmin            ))
       !call check(nf90_put_att(ncid, nf90_global,"YORIG"    , grid%ymin            ))
       !call check(nf90_put_att(ncid, nf90_global,"XCELL"    , grid%dx              ))
       !call check(nf90_put_att(ncid, nf90_global,"YCELL"    , grid%dy              ))
       !call check(nf90_put_att(ncid, nf90_global,"VGTYP"    , -9999                ))
       !call check(nf90_put_att(ncid, nf90_global,"VGTOP"    , 5000.                ))
       !call check(nf90_put_att(ncid, nf90_global,"VGLVLS"   , [1., 0.9938147 ]     ))
       !call check(nf90_put_att(ncid, nf90_global,"GDNAM"    , grid%gName           ))
       !call check(nf90_put_att(ncid, nf90_global,"UPNAM"    , "OUTCM3IO"           ))
       !call check(nf90_put_att_any(ncid, nf90_global,"VAR-LIST",nf90_char, nvars*16, adjustl(var_list_string)))
       !call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "Fire emission file" ))
       !call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
  
     call check(nf90_close(ncid))
     !Cierro NetCDF outFile

     deallocate(var_list)       !Libero memoria
     deallocate(var_desc)       !Libero memoria
     deallocate(file_list)      !Libero memoria
     deallocate(CTS)            !Libero memoria

end subroutine build_CT3

subroutine build_LAI(g,p,laivfile)
   implicit none
   type(grid_type) ,intent(in) :: g
   type(proj_type) ,intent(in) :: p
   character(len=200) :: laivfile
   real, allocatable :: LAIv(:,:,:)
   integer :: date_time_dim_id,var_dim_id,lev_id,tstep_dim_id,var_id,ncid
   character(len=6) :: outfile
   character(len=25),allocatable :: var_list(:),var_desc(:)
   character(len=200),allocatable :: file_list(:)
   integer :: nvars
   character(len=2)::kk

   outfile='LAI3.nc'
   nvars=12

   allocate(LAIv(grid%nx,grid%ny,nvars))  

   !Levanto netcdf input files
   do k=1,nvars
        print*,"k=",k
       call check(nf90_open( trim(file_list(k)), nf90_write, ncid   ))
       call check(   nf90_inq_varid(ncid,'Band1', var_id     ))
       call check(   nf90_get_var(ncid,  var_id , CTS(:,:,k)  ))
       call check(nf90_close(ncid                             ))
   enddo

   var_list =(/ "shrub_frac             ","crop_frac              ","grass_frac             ","nl_tree_frac           ","trop_tree_frac         ","bl_tree_frac           "/)
   var_desc =(/ "shrub fraction         ","crop fraction          ","grass fraction         ","needle tree fraction   ","tropical tree fraction ","broadleaf tree fraction"/) 
   ! Create the NetCDF file
   call check(nf90_create(outFile, NF90_CLOBBER, ncid))
       ! Defino dimensiones
       call check(nf90_def_dim(ncid, "TSTEP"    ,   1    , tstep_dim_id    ))
       call check(nf90_def_dim(ncid, "DATE_TIME",   2    , date_time_dim_id))
       call check(nf90_def_dim(ncid, "COL"      , grid%nx, col_dim_id      ))
       call check(nf90_def_dim(ncid, "ROW"      , grid%ny, row_dim_id      ))
       call check(nf90_def_dim(ncid, "LAY"      ,   1    , lay_dim_id      ))
       call check(nf90_def_dim(ncid, "VAR"      , nvars  , var_dim_id      ))

       !Defino variables
       call check(nf90_def_var(ncid,"TFLAG",NF90_FLOAT    , [date_time_dim_id,var_dim_id,tstep_dim_id], var_id))
       call check(nf90_put_att(ncid, var_id, "units"      , "<YYYYDDD,HHMMSS>" ))
       call check(nf90_put_att(ncid, var_id, "long_name"  , "TFLAG           " ))
       call check(nf90_put_att(ncid, var_id, "var_desc"   , "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                "))

       do k=1, nvars
       write(kk,'I02'),k

         call check(nf90_def_var(ncid, "LAI"//kk ,NF90_FLOAT , [col_dim_id,row_dim_id,lay_dim_id,tstep_dim_id], var_id)) !
         call check(nf90_put_att(ncid, var_id,"units"    , "nondimension    " ))
         call check(nf90_put_att(ncid, var_id,"long_name", trim(var_list(k))  ))
         call check(nf90_put_att(ncid, var_id,"var_desc" , trim(var_desc(k))  ))
       end do
    call check(nf90_enddef(ncid))
    !End NetCDF define mode

    !Abro NetCDF outFile
    call check(nf90_open(outFile, nf90_write, ncid       ))
     do k=1, nvars       
     write(kk,'I0.2'),k
       call check(nf90_inq_varid(ncid,"LAI"//kk ,var_id))               
       call check(nf90_put_var(ncid, var_id, CTS(:,:,6)*(1-CTS(:,:,4)) *(1-CTS(:,:,5)) ))
       !TFLAG:
       call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
       call check(nf90_put_var(ncid, var_id, (/0000000,000000 /) ))





end subroutine





end program


!ESTRUCTURA GENERAL de MEGAN_PREP_Code_Jan_2022:
!(1) Leen Namelist
!(2) call wrf_file         (leen wrfout para sacar algunos parametros de grilla y proyeccion)
!(3) call  megan2_bioemiss (creo que aca es donde interpola)
!(4) call write_*
!       prepmegan4cmaq_cantype.f90:      write(99,'(a)')"CID,ICELL,JCELl,NEEDL,TROPI,BROAD,SHRUB,HERB,CROP"
!       prepmegan4cmaq_grwform.bck.f90:  write(99,'(a)')"gridID,TreeFrac,CropFrac,ShrubFrac,HerbFrac"
!       prepmegan4cmaq_grwform.f90:      write(99,'(a)')"gridID,TreeFrac,CropFrac,ShrubFrac,HerbFrac"

!       prepmegan4cmaq_lai.f90:          write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,LAI01,LAI02,LAI03,LAI04,LAI05,LAI06,LAI07,LAI08,LAI09,LAI10, &

!       prepmegan4cmaq_ecotype.f90:      write(99,'(a)')"gridID,EcotypeID,EcotypeFrac"

!
!       prepmegan4cmaq_ef.f90:           write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,DSRAD,DTEMP,ISOP,MYRC,SABI,LIMO,A_3CAR,OCIM,BPIN,APIN,OMTP,FARN,BCAR,OSQT,MBO,MEOH,ACTO,CO,NO,BIDER,STRESS,OTHER"
!       prepmegan4cmaq_fert.f90:         write(99,888)"CELL_ID,X,Y,LAT,LONG",fertday
!       prepmegan4cmaq_landtype.f90:     write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,LANDTYP"
!       

!       prepmegan4cmaq_nitrogen.f90:     write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,NITROGEN01,NITROGEN02,NITROGEN03,NITROGEN04,NITROGEN05,NITROGEN06,NITROGEN07,NITROGEN08,NITROGEN09,NITROGEN10, &
!       prepmegan4cmaq_arid.f90:         write(99,'(a)')"CID, ICELL, JCELL,  ARID"
!       prepmegan4cmaq_non_arid.f90:     write(99,'(a)')"CID, ICELL, JCELL,  NON_ARID"
!       prepmegan4cmaq_pft.f90:          write(99,'(a)')"CID, ICELL, JCELL,  NT_EG_TEMP,  NT_DC_BORL,  NT_EG_BORL, &
!       prepmegan4cmaq_w126.f90:         write(99,'(a)')"CID, ICELL, JCELL,  W126(ppm-hours)"

!(5) clean-up (deallocate)
