!---------------------------------------------------------------
! author: Ramiro A. Espada. April 2023.
! Based on Prep_code &  MEGEFP32 (UCI-BAI-MEGAN)
!---------------------------------------------------------------
program prepmegan4cmaq

  use netcdf
  implicit none

  type proj_type
     character(16)    :: pName     ! nombre de la proyeccion
     character(7)     :: typ_str   ! String  code for projection TYPE
     integer          :: typ       ! Integer code for projection TYPE
     real             :: alp,bet,gam,xcent,ycent!ref_lat,ref_lon,truelat1,truelat2,stand_lon,pole_lat,pole_lon
     character(125)   :: proj4     ! PROJ4 srs definition.
  end type proj_type

  type grid_type
      character(12)    :: gName     !grid-name
      integer          :: nx,ny,nz  !number of cells in x-y direction (ncols, nrows, nlevs)
      real             :: dx,dy     !x-y cell dimension (x_cell, y_cell)
      real             :: xmin,ymin,xmax,ymax,xc,yc
      real             :: lonmin,latmin,lonmax,latmax
   end type grid_type

  type(proj_type) :: proj
  type(grid_type) :: grid

  integer :: status,iostat
  integer :: i,j,k

  character(len=17):: start_date,end_date
  character(200)   :: griddesc_file,gridname,crop_frac_file,grass_frac_file,shrub_frac_file,tree_frac_file,nl_tree_frac_file,bl_tree_frac_file,tp_tree_frac_file,ecotype_file,laiv_files,GtEcoEF_file,arid_file, narid_file, lt_file,nitro_files

  namelist/control/start_date,end_date,griddesc_file,gridname,crop_frac_file,grass_frac_file,shrub_frac_file,tree_frac_file,nl_tree_frac_file, bl_tree_frac_file,tp_tree_frac_file,ecotype_file,laiv_files,GtEcoEF_file,arid_file, narid_file, lt_file,nitro_files

  !Leo namelist:
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'prepmegan4cmaq: failed to read namelist; error = ',iostat
    stop
  end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file, proj, grid)    !(!) mejorar esta funcion basado en lo que haga IOAPI

  !`MEGAN_CTS` (*Canopy Type Fractions*) 
  call build_CT3(grid,proj,crop_frac_file,tree_frac_file,grass_frac_file,shrub_frac_file,nl_tree_frac_file,tp_tree_frac_file)
  
  !`MEGAN_LAI` (Leaf Area Index).
  !call build_LAIv(grid,proj,laiv_files)

  !`MEGAN_EFS` (emission factors) & `MEGAN_LDF` (*Light Dependence Fractions*) 
  call build_EFS_LDF(grid,proj,GtEcoEF_file,ecotype_file,crop_frac_file,tree_frac_file,grass_frac_file,shrub_frac_file)
  
  !!Archivos para BDSNP:
  !call BDSNP_AFILE()   !int Arid (0/1) &  call BDSNP_NAFILE()  !int Non-Arid (0/1) &  call BDSNP_LFILE()   !int Land types (1:24)
  call BDSNP_LAND(grid,proj,arid_file,narid_file,lt_file)

  call BDSNP_NFILE(grid,proj,nitro_files)   !float nitrogeno01, nitrogeno02,...,nitrogeno12  monthly nitrogen deposition in ng of N /m2/s
  !call BDSNP_FFILE()                       !float fert01,fert02,...,fert  daily fertilizer aplication. unit: ng of N/m2 

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

! !Interfaz a "date"
! function date(date_str, fmt_str) result(output)
!   implicit none
!   character(*), intent(in) :: date_str, fmt_str
!   character(256)           :: command
!   character(20)            :: output
!   command="date -d "//trim(date_str)//" '+"//trim(fmt_str)//"'  > tmp_date.txt"
!   call system( trim(command) )
!   !print*,trim(command)
!   open(9, file='tmp_date.txt', status='old',action='read'); read(9, '(A)', iostat=status) output;  close(9)
!   call system('rm tmp_date.txt')
! end function

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
        open(2,file=griddescFile, status='old', action='read')                   !GRIDDESC:
           read(2,*) dummyvar;                                                   !' '
           read(2,*) p%pName;                                                    !projName
           read(2,*) p%typ,p%alp,p%bet,p%gam,p%xcent,p%ycent                     !map_proj truelat1 truelat2 stand_lon ref_lon ref_lat
           read(2,*) dummyvar;                                                   !' '
           read(2,*) g%gName;                                                    !gridName
           read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%nx,g%ny                   !projName xorig yorig xcell ycell nrows ncols
        close(2)
 end subroutine

 subroutine createNetCDF(outFile,p,g,var_list,var_unit,var_desc)
    implicit none
    type(grid_type) , intent(in) :: g
    type(proj_type) , intent(in) :: p
    character(len=10), intent(in) :: outFile
    character(len=25), allocatable:: var_list(:),var_unit(:),var_desc(:)
    integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,var_id
    integer :: nvars
    character(800) :: var_list_string
 
    nvars=size(var_list)
    write(var_list_string,*) var_list                !este es un global attr importante.
    
    call check(nf90_create(outFile, NF90_CLOBBER, ncid))
        ! Defino dimensiones
        call check(nf90_def_dim(ncid, "TSTEP"    ,   1    , tstep_dim_id     ))
        call check(nf90_def_dim(ncid, "DATE_TIME",   2    , date_time_dim_id ))
        call check(nf90_def_dim(ncid, "COL"      , g%nx   , col_dim_id       ))
        call check(nf90_def_dim(ncid, "ROW"      , g%ny   , row_dim_id       ))
        call check(nf90_def_dim(ncid, "LAY"      ,   1    , lay_dim_id       ))
        call check(nf90_def_dim(ncid, "VAR"      , nvars  , var_dim_id       ))
        !Defino variables
        call check(nf90_def_var(ncid,"TFLAG",NF90_FLOAT    , [date_time_dim_id,var_dim_id,tstep_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id, "units"      , "<YYYYDDD,HHMMSS>" ))
        call check(nf90_put_att(ncid, var_id, "long_name"  , "TFLAG           " ))
        call check(nf90_put_att(ncid, var_id, "var_desc"   , "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                "))
        do k=1, nvars
          call check(nf90_def_var(ncid, trim(var_list(k)) , NF90_FLOAT, [col_dim_id,row_dim_id,lay_dim_id,tstep_dim_id], var_id)) 
          call check(nf90_put_att(ncid, var_id,"units"    , trim(var_unit(k)) ))
          call check(nf90_put_att(ncid, var_id,"long_name", trim(var_list(k)) ))
          call check(nf90_put_att(ncid, var_id,"var_desc" , trim(var_desc(k)) ))
        end do
        ! Defino attributos
        call check(nf90_put_att(ncid, nf90_global,"IOAPI_VERSION", "ioapi-3.2: \$Id: init3" ))
        call check(nf90_put_att(ncid, nf90_global,"EXEC_ID"  , "????????????????"   ))
        call check(nf90_put_att(ncid, nf90_global,"FTYPE"    , 1                    ))
        call check(nf90_put_att(ncid, nf90_global,"SDATE"    , 2023001              ))!stat_date (int)
        call check(nf90_put_att(ncid, nf90_global,"STIME"    , 000000               ))
        call check(nf90_put_att(ncid, nf90_global,"WDATE"    , 2023001              ))
        call check(nf90_put_att(ncid, nf90_global,"WTIME"    , 000000               ))
        call check(nf90_put_att(ncid, nf90_global,"CDATE"    , 2023001              ))
        call check(nf90_put_att(ncid, nf90_global,"CTIME"    , 000000               ))
        call check(nf90_put_att(ncid, nf90_global,"TSTEP"    , 10000                ))
        call check(nf90_put_att(ncid, nf90_global,"NTHIK"    , 1                    ))!no sé que es.
        call check(nf90_put_att(ncid, nf90_global,"NCOLS"    , g%nx                 ))
        call check(nf90_put_att(ncid, nf90_global,"NROWS"    , g%ny                 ))
        call check(nf90_put_att(ncid, nf90_global,"NLAYS"    , 1                    ))!grid%nz
        call check(nf90_put_att(ncid, nf90_global,"NVARS"    , nvars                ))
        call check(nf90_put_att(ncid, nf90_global,"GDTYP"    , p%typ                ))
        call check(nf90_put_att(ncid, nf90_global,"P_ALP"    , p%alp                ))
        call check(nf90_put_att(ncid, nf90_global,"P_BET"    , p%bet                ))
        call check(nf90_put_att(ncid, nf90_global,"P_GAM"    , p%gam                ))
        call check(nf90_put_att(ncid, nf90_global,"XCENT"    , p%xcent              ))
        call check(nf90_put_att(ncid, nf90_global,"YCENT"    , p%ycent              ))
        call check(nf90_put_att(ncid, nf90_global,"XORIG"    , g%xmin               ))
        call check(nf90_put_att(ncid, nf90_global,"YORIG"    , g%ymin               ))
        call check(nf90_put_att(ncid, nf90_global,"XCELL"    , g%dx                 ))
        call check(nf90_put_att(ncid, nf90_global,"YCELL"    , g%dy                 ))
        call check(nf90_put_att(ncid, nf90_global,"VGTYP"    , -9999                ))!no sé que es.
        call check(nf90_put_att(ncid, nf90_global,"VGTOP"    , 5000.                ))!no sé que es.
        call check(nf90_put_att(ncid, nf90_global,"VGLVLS"   , [1., 0.9938147 ]     ))!no sé que es.
        call check(nf90_put_att(ncid, nf90_global,"GDNAM"    , g%gName              ))
        call check(nf90_put_att(ncid, nf90_global,"UPNAM"    , "OUTCM3IO"           ))!no sé que es.
        call check(nf90_put_att_any(ncid, nf90_global,"VAR-LIST",nf90_char, nvars*16, adjustl(var_list_string)))
        call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "MEGAN input file"   ))
        call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))

     call check(nf90_enddef(ncid))
     !End NetCDF define mode
 end subroutine createNetCDF
 
 function get2DvarFromNetCDF(ncfile, varname, nx, ny)  result(GRIDVAR)
        implicit none
        character(50), intent(in) :: ncfile
        character(5) , intent(in) :: varname
        integer, intent(in)       :: nx,ny
        real, allocatable :: GRIDVAR(:,:)
        integer var_id,ncid
 
        allocate(GRIDVAR(nx,ny))

        call check(nf90_open(trim(ncfile), nf90_write, ncid ))
             call check(   nf90_inq_varid(ncid,trim(varname), var_id   ))
             call check(   nf90_get_var(ncid, var_id , GRIDVAR   ))
        call check(nf90_close(ncid))
 end function get2DvarFromNetCDF

 function get2DvarFromNetCDFint(ncfile, varname, nx, ny)  result(GRIDVAR)
        implicit none
        character(50), intent(in) :: ncfile
        character(5) , intent(in) :: varname
        integer, intent(in)       :: nx,ny
        integer, allocatable :: GRIDVAR(:,:)
        integer var_id,ncid
 
        allocate(GRIDVAR(nx,ny))

        call check(nf90_open(trim(ncfile), nf90_write, ncid ))
             call check(   nf90_inq_varid(ncid,trim(varname), var_id   ))
             call check(   nf90_get_var(ncid, var_id , GRIDVAR   ))
        call check(nf90_close(ncid))
 end function get2DvarFromNetCDFint
 !----------------------------------
 !  MEGAN_CTS 
 !---------------------------------
 subroutine build_CT3(g,p,cropfile,treefile,grassfile,shrubfile,nltreefile,troptreefile)
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(len=200) :: cropfile,treefile,grassfile,shrubfile,nltreefile,troptreefile
    real, allocatable :: CTS(:,:,:)
    character(len=10) :: outfile
    character(len=25),allocatable :: var_list(:),var_desc(:),var_unit(:)
    integer :: nvars
    integer :: var_id,ncid

 
    print*,"Building MEGAN_CTS file ..."
 
    outfile='CT3.nc'
    nvars=6
 
    allocate(var_list(nvars))  
    allocate(var_unit(nvars))  
    allocate(var_desc(nvars))  
    allocate(CTS(g%nx,g%ny,nvars))  
 
    CTS(:,:,1)=get2DvarFromNetCDF(    treefile, "Band1", g%nx, g%ny)
    CTS(:,:,2)=get2DvarFromNetCDF(   shrubfile, "Band1", g%nx, g%ny)
    CTS(:,:,3)=get2DvarFromNetCDF(    cropfile, "Band1", g%nx, g%ny)
    CTS(:,:,4)=get2DvarFromNetCDF(   grassfile, "Band1", g%nx, g%ny)
    CTS(:,:,5)=get2DvarFromNetCDF(  nltreefile, "Band1", g%nx, g%ny)
    CTS(:,:,6)=get2DvarFromNetCDF(troptreefile, "Band1", g%nx, g%ny)
    where ( CTS< 0.0 )
         CTS=0.0
    end where
 
    var_list=(/ 'NEEDL      '     ,'TROPI      ' ,'BROAD      ','SHRUB      ','HERB       ','CROP       '/)
    var_desc=(/ "needle tree fraction   ","tropical tree fraction ","broadleaf tree fraction", "shrub fraction         ","crop fraction          ","grass fraction         " /) 
    var_unit=spread("nondimension",1,nvars)
 
    ! Create the NetCDF file
    call createNetCDF(outFile,p,g,var_list,var_unit,var_desc)
 
    !Abro NetCDF outFile
     call check(nf90_open(outFile, nf90_write, ncid       ))
        !shrub
        call check(nf90_inq_varid(ncid,"SHRUB" ,var_id))
        call check(nf90_put_var(ncid, var_id, CTS(:,:,2)  ))        
        !crop
        call check(nf90_inq_varid(ncid,"CROP"  ,var_id))
        call check(nf90_put_var(ncid, var_id, CTS(:,:,3)  ))        
        !grass
        call check(nf90_inq_varid(ncid,"HERB"  ,var_id))
        call check(nf90_put_var(ncid, var_id, CTS(:,:,4)  ))        
        !needleleaf tree
        call check(nf90_inq_varid(ncid,"NEEDL" ,var_id))
        call check(nf90_put_var(ncid, var_id, CTS(:,:,1) * (1.0-CTS(:,:,6)) * CTS(:,:,5)     )) 
        !boradleaf tree
        call check(nf90_inq_varid(ncid,"BROAD" ,var_id))               
        call check(nf90_put_var(ncid, var_id, CTS(:,:,1) * (1.0-CTS(:,:,6)) * (1.0-CTS(:,:,5)) ))
        !tropical tree
        call check(nf90_inq_varid(ncid,"TROPI" ,var_id))
        call check(nf90_put_var(ncid, var_id, CTS(:,:,1)*CTS(:,:,6)                    )) 
        !TFLAG:
        call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
        call check(nf90_put_var(ncid, var_id, (/0000000,000000 /) ))
      call check(nf90_close(ncid))
      !Cierro NetCDF outFile
 
      deallocate(CTS)            !Libero memoria
 
 end subroutine build_CT3
 
 !----------------------------------
 !  MEGAN_LAI 
 !---------------------------------
 subroutine build_LAIv(g,p,laivfile)
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(len=200) :: laivfile
    real, allocatable :: LAIv(:,:,:)
    character(len=25),allocatable :: var_list(:),var_desc(:),var_unit(:)
    character(len=10) :: outfile
    integer ::var_id,ncid
    integer :: nvars
    character(len=2)::kk
    
    print*,"Building MEGAN_LAI file ..."
    
    outfile='LAI3.nc'
    nvars=12
 
    allocate(var_list(nvars))  
    allocate(var_unit(nvars))  
    allocate(var_desc(nvars))  
    allocate(LAIv(g%nx,g%ny,nvars))  
 
    !Levanto netcdf input files
    do k=1,nvars
        write(kk,'(I0.2)') k
        LAIv(:,:,k)=get2DvarFromNetCDF(trim(laivfile)//kk//".nc", "Band1", g%nx, g%ny)
    enddo
    where (LAIv < 0.0 )
            LAIv=0.0
    endwhere

    var_list=(/"LAI01","LAI02","LAI03","LAI04","LAI05","LAI06","LAI07","LAI08","LAI09","LAI10","LAI11","LAI12" /)
    var_desc=(/"LAI01","LAI02","LAI03","LAI04","LAI05","LAI06","LAI07","LAI08","LAI09","LAI10","LAI11","LAI12" /)
    var_unit=spread("nondimension",1,nvars)
 
    !Creo NetCDF file
    call createNetCDF(outFile,p,g,var_list,var_unit,var_desc)
    !Abro NetCDF outFile
    call check(nf90_open(outFile, nf90_write, ncid       ))
      do k=1, nvars       
        write(kk,'(I0.2)') k
        call check(nf90_inq_varid(ncid,"LAI"//kk ,var_id))               
        call check(nf90_put_var(ncid, var_id, LAIv(:,:,k)/1000.0 ))
      enddo
        !TFLAG:
        call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
        call check(nf90_put_var(ncid, var_id, (/0000000,000000 /) ))
     !Cierro NetCDF outFile
     call check(nf90_close( ncid ))
 
 end subroutine
 
  !----------------------------------
  !  MEGAN_LDF & MEGAN_EF
  !---------------------------------
   subroutine build_EFS_LDF(g,p,GtEcoEF_file,ecotypefile,cropfile,treefile,grassfile,shrubfile)
      implicit none
      type(grid_type) ,intent(in) :: g
      type(proj_type) ,intent(in) :: p
      character(len=50),intent(in) :: ecotypefile,cropfile,treefile,grassfile,shrubfile,GtEcoEF_file
      character(len=10) :: outfileEF,outfileLDF
      integer i,j,k
      integer :: var_id,ncid
      
      real,    allocatable :: GTYP(:,:,:)       !growth type fraction
      integer, allocatable :: ECOTYPE(:,:)      !este puede ser int tamb
      integer :: EcoID
      real, allocatable :: OUTGRID(:,:,:)

      character(len=5) :: GTYP_LIST(4) ! crop, tree, grass, shrub
      character(len=5) :: GtID
      character(len=25), allocatable :: var_list(:),var_desc(:),var_unit(:) !19 EF + 4 LDF
      integer :: nvars
      real  :: EF(23)

      print*,"Building MEGAN_EFS & MEGAN_LDF ..."
      
      nvars=23
       outfileEF='EFMAP.nc'
      outfileLDF='LDF.nc'
 
      allocate( ECOTYPE(g%nx, g%ny   ))   !ecotype         
      allocate(    GTYP(g%nx, g%ny,4 ))   !growthtype fracs

      !Levanto netcdf input files
      call check(nf90_open(trim(ecotypefile), nf90_write, ncid ))
           call check(   nf90_inq_varid(ncid,'Band1', var_id   ))
           call check(   nf90_get_var(ncid, var_id , ECOTYPE   ))
      call check(nf90_close(ncid ))

      GTYP_LIST=(/'crop ','tree ','herb ','shrub'/)

      GTYP(:,:,1)=get2DvarFromNetCDF( cropfile, "Band1", g%nx, g%ny)
      GTYP(:,:,2)=get2DvarFromNetCDF( treefile, "Band1", g%nx, g%ny)
      GTYP(:,:,3)=get2DvarFromNetCDF(grassfile, "Band1", g%nx, g%ny)
      GTYP(:,:,4)=get2DvarFromNetCDF(shrubfile, "Band1", g%nx, g%ny)
      where (GTYP < 0.0 )
              GTYP=0.0
      endwhere
      GTYP=GTYP/100  !las frac están en porcentaje.

     !Arranco con una tabla: GtEcoEF_file:
     !growtypeId   ecotypeID   var1,   var2, ... ,  var19,  var20, ... ,  var23
     !crop         1           EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
     !crop         2           EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
     !crop         3           EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
     !crop         4           EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
     !....         ...         EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
     !tree         1700        EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
 
     allocate(OUTGRID(g%nx, g%ny, 23))   !outgrids EF1,EF2,...,LDF1,LDF2,..
     OUTGRID=0.0
     j=0
     open(unit=1,file=trim(GtEcoEF_file),status='old',action='read')
       iostat=0
       do while(iostat == 0)       !loop por cada fila
          read(1,*,iostat=iostat)GtID,EcoID,EF
          
          if (j /= FINDLOC(GTYP_LIST, GtID,1) ) then
                j=FINDLOC(GTYP_LIST, GtID,1)
                print*,"   Processing Growth-type: "//GtID
          endif
          !=======> (!) ACÁ está el cuello de botella <=====
          do i=1,23    !nvars: EF/LDF
            where ( ECOTYPE == EcoID )
                OUTGRID(:,:,i) = OUTGRID(:,:,i) + GTYP(:,:,j) * EF(i)
            endwhere
          enddo !i (var)
          !=======> (!) ACÁ está el cuello de botella <=====
       enddo !each row.
     close(1)
 
     deallocate(GTYP)
     deallocate(ECOTYPE)
      
     allocate(var_list(23),var_desc(23),var_unit(23))

     var_list=(/"EF_ISOP   ","EF_MBO    ","EF_MT_PINE","EF_MT_ACYC","EF_MT_CAMP","EF_MT_SABI","EF_MT_AROM","EF_NO     ","EF_SQT_HR ","EF_SQT_LR ", "EF_MEOH   ","EF_ACTO   ","EF_ETOH   ","EF_ACID   ","EF_LVOC   ","EF_OXPROD ","EF_STRESS ","EF_OTHER  ","EF_CO     ", "LDF03     ", "LDF04     ", "LDF05     ", "LDF06     " /)
     var_desc=var_list
     var_unit=spread("nanomol/m^2/s",1,nvars)
   
     !EF.nc ==============================
     ! Create the NetCDF file
     call createNetCDF(outFileEF,p,g,var_list,var_unit,var_desc)
     !Abro NetCDF outFile
     call check(nf90_open(outFileEF, nf90_write, ncid       ))
      do k=1, 23       
        call check(nf90_inq_varid(ncid,var_list(k) ,var_id ))
        call check(nf90_put_var(ncid, var_id, OUTGRID(:,:,k)))
      enddo
     call check(nf90_close( ncid ))
     !Cierro NetCDF outFile================
     !!LDF.nc ==============================
     !! Create the NetCDF file
     !call createNetCDF(outFileLDF,p,g,var_list,var_unit,var_desc)
     !!Abro NetCDF outFile
     !call check(nf90_open(outFileLDF, nf90_write, ncid       ))
     ! do k=20,23       
     !   call check(nf90_inq_varid(ncid,var_list(k) ,var_id))               
     !   call check(nf90_put_var(ncid, var_id, OUTGRID(:,:,k) ))
     ! enddo
     !call check(nf90_close( ncid ))
     !!Cierro NetCDF outFile================

     deallocate(OUTGRID)
   end subroutine


 !----------------------------------
 !  BDSNP_ARID, BDSNP_NONARID & BDSNP_LANDTYPE
 !---------------------------------
  subroutine BDSNP_LAND(g,p,arid_file, narid_file,lt_file)
     implicit none
     type(grid_type) ,intent(in) :: g
     type(proj_type) ,intent(in) :: p
     character(len=200),intent(in) :: arid_file,narid_file,lt_file
     real, allocatable :: LANDGRID(:,:,:)
     character(len=25),allocatable :: var_list(:),var_desc(:),var_unit(:)
     integer :: var_id,ncid
     integer :: nvars
     integer :: lt
     character(len=2) ::landtypeId
     character(len=10) :: outfile

     print*,"Building BDSNP_ARID, BDSNP_NONARID & BDSNP_LANDTYPE ..."
     outfile='LAND.nc'!(/"ARID.nc    ", "LANDTYPE.nc", "NONARID.nc "/) !
     nvars=3
     allocate(LANDGRID(g%nx,g%ny,nvars))
     allocate(var_list(nvars))
     allocate(var_unit(nvars))
     allocate(var_desc(nvars))

     LANDGRID(:,:,1) = get2DvarFromNetCDFint( arid_file, 'Band1', g%nx, g%ny)
     LANDGRID(:,:,2) = get2DvarFromNetCDFint(narid_file, 'Band1', g%nx, g%ny)
     
     LANDGRID(:,:,3) =0.0 
     do lt=1,24
             write(landtypeId, '(I0.2)') lt
             LANDGRID(:,:,3) = LANDGRID(:,:,3)+lt* get2DvarFromNetCDFint(trim(lt_file)//landtypeId//".nc", 'Band1', g%nx, g%ny)
     enddo

     var_list=(/ 'ARID    ', 'NONARID ', 'LANDTYPE' /)
     var_unit=spread('1 or 0      ',1,nvars) 
     var_desc=(/ 'ARID    ', 'NONARID ', 'LANDTYPE' /)
      
     ! Create the NetCDF file
     call createNetCDF(outFile,p,g,var_list,var_unit,var_desc)

     !Abro NetCDF outFile
     call check(nf90_open(outFile, nf90_write, ncid       ))
     do k=1, nvars    
        call check(nf90_inq_varid(ncid,var_list(k) ,var_id ))
        call check(nf90_put_var(ncid, var_id, LANDGRID(:,:,k)))
     enddo
     call check(nf90_close( ncid ))
  end subroutine BDSNP_LAND

  subroutine BDSNP_NFILE(g,p,nitro_files)   
     implicit none
     type(grid_type) ,intent(in) :: g
     type(proj_type) ,intent(in) :: p
     character(len=200),intent(in) :: nitro_files
     real, allocatable :: NITRO(:,:,:)
     character(len=10) :: outfile
     character(len=25),allocatable :: var_list(:),var_desc(:),var_unit(:)
     integer :: var_id,ncid
     integer :: nvars
     integer :: k
     character(len=2) ::kk
                                                                            
     print*,"Building BDSNP_NFILE ..."
     outfile='NDEP.nc'
     nvars=12

     allocate(var_list(nvars))
     allocate(var_unit(nvars))
     allocate(var_desc(nvars))
     allocate(NITRO(g%nx,g%ny,nvars))  
                                                                             
     !Levanto netcdf input files
     do k=1,nvars
         write(kk,'(I0.2)') k
         call check(nf90_open(trim(nitro_files)//kk//".nc", nf90_write, ncid ))
              call check(   nf90_inq_varid(ncid,'Band1', var_id ))
              call check(   nf90_get_var(ncid, var_id , NITRO(:,:,k)  ))
              where (NITRO < 0.0 )
                      NITRO=0.0
              endwhere
         call check(nf90_close(ncid ))
     enddo
     NITRO=NITRO*1E+12! convert from kg/m2/s to ng/m2/s

     var_list=(/ 'NITROGEN01','NITROGEN02','NITROGEN03','NITROGEN04','NITROGEN05','NITROGEN06','NITROGEN07','NITROGEN08','NITROGEN09','NITROGEN10','NITROGEN11','NITROGEN12' /)
    var_unit = spread( "ng/m2/s         ",1,nvars)
    var_desc=spread("monthly average total nitrogen deposition",1,nvars) 
                                                                
     ! Create the NetCDF file
     call createNetCDF(outFile,p,g,var_list,var_unit,var_desc)
                                                                
     !Abro NetCDF outFile
     call check(nf90_open(outFile, nf90_write, ncid       ))
      do k=1, nvars    
        call check(nf90_inq_varid(ncid,var_list(k) ,var_id ))
        call check(nf90_put_var(ncid, var_id, NITRO(:,:,k)))
      enddo
     call check(nf90_close( ncid ))

  end subroutine 
end program
