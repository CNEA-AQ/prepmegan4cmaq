!---------------------------------------------------------------
! author: Ramiro A. Espada. April 2023.
! Based on UCI-BAI-MEGAN Prep_code and megan_bio_emis.f90 (NCAR)
!---------------------------------------------------------------
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
  character(200) :: griddesc_file, cropf_file,grassf_file,shrubf_file,treef_file,nl_treef_file,bl_treef_file,tp_treef_file,ecotype_file,laiv_files

  namelist/control/start_date,end_date,griddesc_file, cropf_file,grassf_file,shrubf_file,treef_file,nl_treef_file, bl_treef_file,tp_treef_file,ecotype_file,laiv_files!,lai_file,veg_cov_file,wrfout_file

  !Leo namelist:
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'prepmegan4cmaq: failed to read namelist; error = ',iostat
    stop
  end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file, proj, grid)    !(!) TO-DO: mejorar esta funcion basado en lo que haga IOAPI

  !`MEGAN_CTS` (*Canopy Type Fractions*) 
  call build_CT3(grid,proj,cropf_file,treef_file,grassf_file,shrubf_file,nl_treef_file,tp_treef_file)
  
  !`MEGAN_LAI` (Leaf Area Index).
  !call build_LAIv(grid,proj,laiv_files)

  !!`MEGAN_LDF` (*Light Dependence Fractions*)
  !call build_LDF()

  !!`MEGAN_EFS` (emission factors).
  !call build_EFS()

  !!Archivos para BDSNP: (Solo es regrillar netcdf):
  !call BDSNP_AFILE()   !int Arid (0/1) , int Non-Arid (0/1)

  !call BDSNP_LFILE()   !int Land types (1:24)

  !call BDSNP_FFILE()   !float fert01,fert02,...,fert  daily fertilizer aplication. unit: ng of N/m2 

  !call BDSNP_NFILE()   !float nitrogeno01, nitrogeno02,...,nitrogeno12  monthly nitrogen deposition in ng of N /m2/s

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
        open(2,file=griddescFile, status='old', action='read')                                  !GRIDDESC:
           read(2,*) dummyvar;                                                   !' '
           read(2,*) p%pName;                                                    !projName
           read(2,*) p%typ,p%truelat1,p%truelat2,p%stand_lon,p%ref_lon,p%ref_lat !map_proj truelat1 truelat2 stand_lon ref_lon ref_lat
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
    integer :: date_time_dim_id,var_dim_id,lev_id,tstep_dim_id,var_id,ncid
    integer :: nvars
    real :: missing_value
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
        call check(nf90_put_att(ncid, nf90_global,"SDATE"    , 2023001              ))   !int
        call check(nf90_put_att(ncid, nf90_global,"STIME"    , 000000               ))
        call check(nf90_put_att(ncid, nf90_global,"WDATE"    , 2023001              ))
        call check(nf90_put_att(ncid, nf90_global,"WTIME"    , 000000               ))
        call check(nf90_put_att(ncid, nf90_global,"CDATE"    , 2023001              ))
        call check(nf90_put_att(ncid, nf90_global,"CTIME"    , 000000               ))
        call check(nf90_put_att(ncid, nf90_global,"TSTEP"    , 10000                ))
        call check(nf90_put_att(ncid, nf90_global,"NTHIK"    , 1                    ))
        call check(nf90_put_att(ncid, nf90_global,"NCOLS"    , g%nx                 ))
        call check(nf90_put_att(ncid, nf90_global,"NROWS"    , g%ny                 ))
        call check(nf90_put_att(ncid, nf90_global,"NLAYS"    , 1                    ))!grid%nz
        call check(nf90_put_att(ncid, nf90_global,"NVARS"    , nvars                ))
        call check(nf90_put_att(ncid, nf90_global,"GDTYP"    , 1                    ))
        call check(nf90_put_att(ncid, nf90_global,"P_ALP"    , -50.                 ))
        call check(nf90_put_att(ncid, nf90_global,"P_BET"    , -20.                 ))
        call check(nf90_put_att(ncid, nf90_global,"P_GAM"    , -65.                 ))
        call check(nf90_put_att(ncid, nf90_global,"XCENT"    , p%ref_lon            ))
        call check(nf90_put_att(ncid, nf90_global,"YCENT"    , p%ref_lat            ))
        call check(nf90_put_att(ncid, nf90_global,"XORIG"    , g%xmin               ))
        call check(nf90_put_att(ncid, nf90_global,"YORIG"    , g%ymin               ))
        call check(nf90_put_att(ncid, nf90_global,"XCELL"    , g%dx                 ))
        call check(nf90_put_att(ncid, nf90_global,"YCELL"    , g%dy                 ))
        call check(nf90_put_att(ncid, nf90_global,"VGTYP"    , -9999                ))
        call check(nf90_put_att(ncid, nf90_global,"VGTOP"    , 5000.                ))
        call check(nf90_put_att(ncid, nf90_global,"VGLVLS"   , [1., 0.9938147 ]     ))
        call check(nf90_put_att(ncid, nf90_global,"GDNAM"    , g%gName              ))
        call check(nf90_put_att(ncid, nf90_global,"UPNAM"    , "OUTCM3IO"           ))
        call check(nf90_put_att_any(ncid, nf90_global,"VAR-LIST",nf90_char, nvars*16, adjustl(var_list_string)))
        call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "Fire emission file" ))
        call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
 
     call check(nf90_enddef(ncid))
     !End NetCDF define mode
 end subroutine createNetCDF
 
 
 
 
 !----------------------------------
 !  MEGAN_CTS 
 !---------------------------------
 subroutine build_CT3(g,p,cropfile,treefile,grassfile,shrubfile,nltreefile,troptreefile)
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(len=200) :: cropfile,treefile,grassfile,shrubfile,nltreefile,troptreefile
    !real, allocatable :: crop(:,:),shrub(:,:),grass(:,:),tree(:,:),nl_tree(:,:),trop_tree(:,:),bl_tree(:,:),laiv(:,:)
    real, allocatable :: CTS(:,:,:)
    integer :: date_time_dim_id,var_dim_id,lev_id,tstep_dim_id,var_id,ncid
    character(len=10) :: outfile
    character(len=25),allocatable :: var_list(:),var_desc(:),var_unit(:)
    character(len=200),allocatable :: file_list(:)
    integer :: nvars
 
    print*,"Building MEGAN_CTS file ..."
 
    outfile='CT3.nc'
    nvars=6
 
    allocate(var_list(nvars))  
    allocate(var_unit(nvars))  
    allocate(var_desc(nvars))  
    allocate(file_list(nvars))  
    allocate(CTS(grid%nx,grid%ny,nvars))  
 
    file_list=(/treefile, shrubfile,cropfile, grassfile, nltreefile, troptreefile/)      !No cambiar, EL ORDEN IMPORTA (do not touch, the order matters)!
    !Levanto netcdf input files
    do k=1,nvars
        call check(nf90_open( trim(file_list(k)), nf90_write, ncid   ))
        call check(   nf90_inq_varid(ncid,'Band1', var_id     ))
        call check(   nf90_get_var(ncid,  var_id , CTS(:,:,k)  ))
        call check(nf90_close(ncid                             ))
    enddo
    where ( CTS< 0.0 )
         CTS=0.0
    end where
 
    var_list=(/'NEEDL      '     ,'TROPI      ' ,'BROAD      ','SHRUB      ','HERB       ','CROP       '/)
    var_desc=(/ "needle tree fraction   ","tropical tree fraction ","broadleaf tree fraction", "shrub fraction         ","crop fraction          ","grass fraction         " /) 
    var_unit=(/ "nondmension","nondmension","nondmension","nondmension","nondmension","nondmension" /)
 
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
    integer :: date_time_dim_id,var_dim_id,lev_id,tstep_dim_id,var_id,ncid
    character(len=25),allocatable :: var_list(:),var_desc(:),var_unit(:)
    character(len=10) :: outfile
    integer :: nvars
    character(len=2)::kk
    real :: missing_value
    
    print*,"Building MEGAN_LAI file ..."
    
    outfile='LAI3.nc'
    nvars=12
 
    allocate(var_list(nvars))  
    allocate(var_unit(nvars))  
    allocate(var_desc(nvars))  
    allocate(LAIv(grid%nx,grid%ny,nvars))  
 
    !Levanto netcdf input files
    do k=1,nvars
        write(kk,'(I0.2)') k
        call check(nf90_open(trim(laivfile)//kk//".nc", nf90_write, ncid ))
             call check(   nf90_inq_varid(ncid,'Band1', var_id ))
             call check(   nf90_get_var(ncid, var_id , LAIv(:,:,k)  ))
             where (LAIv < 0.0 )
                     LAIv=0.0
             endwhere
        call check(nf90_close(ncid ))
    enddo
  
    var_list=(/"LAI01","LAI02","LAI03","LAI04","LAI05","LAI06","LAI07","LAI08","LAI09","LAI10","LAI11","LAI12" /)
    var_desc=(/"LAI01","LAI02","LAI03","LAI04","LAI05","LAI06","LAI07","LAI08","LAI09","LAI10","LAI11","LAI12" /)
    var_unit=(/"nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension" /)
 
    ! Create the NetCDF file
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
 
     call check(nf90_close( ncid ))
     !Cierro NetCDF outFile
 
 end subroutine
 

  !----------------------------------
  !  MEGAN_LDF & MEGAN_EF
  !---------------------------------

  !===============================================================================================================
  !Las siguientes funciones usan MFEP (Megan Factor Emission Processor): acà voy a ir describiendo que entiendo que hace este
  !procesador:
  
  !(1) Carga "bases de datos"
  
  !==> EFv210806.csv <==  Perfiles de EF y LDF de cada VegID
  !VegID,CommonName,GenusGroup,Family,GrowthForm,Type,EF1,EF2,EF3,EF4,EF5,EF6,EF7,EF8,EF9,EF10,EF11,EF12,EF13,EF14,EF15,EF16,EF17,EF18,EF19,LDF3,LDF4,LDF5,LDF6,References,Comment
  !
  !==> SpeciationCrop210806.csv                    <==Esta setiado all ==1  Cuanto de cada veg hay en cada growth-type y cada ecotype
  !EcoTypeID,VegID,<Crop>SpecFrac
  !==> SpeciationHerb210806.csv                    <==Esta setiado all ==1 
  !EcoTypeID,VegID,>Herb>SpecFrac
  !==> SpeciationShrub210806.csv                   <==Esta setiado all ==1 
  !EcoTypeID,VegID,<Shrub>SpecFrac
  !==> SpeciationTree210725.csv                    <==Este si tiene datos
  !EcoTypeID,VegID,<tree>Specfrac,EcoTypeDescription
  !
  !==> grid_ecotype.tx_12km.csv                    <== fraccion de cada ecotypo prepmegan4cmaq_ecotype.f90
  !gridID,EcotypeID,EcotypeFrac                        
  !==> grid_growth_form.tx_12km.csv                <== fraccion de cada growth  prepmegan4cmaq_grwform.f90         
  !gridID,TreeFrac,CropFrac,ShrubFrac,HerbFrac
  !
  !(2) Hace sql-queries para cruzar info
  ! def build_interm_query(growthform,ef_s=1,ef_e=18,ldf_s=3,ldf_e=6):
  !     """
  !     Function to build the intermediate query
  !     """
  !    query_str = f"CREATE TABLE 'Intermediate{growthform}EcoEF' AS \
  !                  SELECT Speciation{growthform}.EcoTypeID, "
  !    for i in range(ef_s,ef_e+1,1):
  !        q = f"
  !        Sum( [EF{i}] * [{growthform}Specfrac] ) AS {growthform}EcoEF{i}, "
  !        query_str+=q
  !    for j in range(ldf_s,ldf_e+1,1):
  !        q2 = f"Sum([LDF{j}]*[{growthform}Specfrac]) AS {growthform}EcoLDF{j}, "
  !        query_str+=q2
  

  !    query_str = query_str.rstrip(', ')
  !    query_str+=f" FROM Speciation{growthform} \
  !                 INNER JOIN EF ON Speciation{growthform}.VegID = EF.VegID \
  !                 GROUP BY Speciation{growthform}.EcoTypeID;"
  !    return query_str
  ! 
  ! 
  ! def Ecotype_Tree_EF(conn, EFa=1, EFz=18, LDFa=3, LDFz=6):
  !     c = conn.cursor()
  !     query_str = build_interm_query("Tree",EFa, EFz, LDFa, LDFz)
  !     c.execute(query_str)
  !     print("'IntermediateTreeEcoEF' Table Created")
  !...(idem con el resto de los growthtypes)... 
  ! 
  !(3) La cuenta seria: (pensarla)
  !
  !def build_final_query(ef_s=1,ef_e=18,ldf_s=3,ldf_e=6):
  !    """
  !    build the final query to calculate the grid EF
  !    """
  !    query_str = "CREATE TABLE 'OutputGridEF' AS \
  !                SELECT GridGrowthForm.gridID, "
  !    for i in range(ef_s,ef_e+1,1):
  !        q = f"Sum([EcotypeFrac]*(([CropFrac]*[CropEcoEF{i}])\
  !        +([TreeFrac]*[TreeEcoEF{i}])\
  !        +([HerbFrac]*[HerbEcoEF{i}])\
  !        +([ShrubFrac]*[ShrubEcoEF{i}]))) AS EF{i}, "
  !        query_str += q
  !
  !    for j in range(ldf_s,ldf_e+1,1):
  !        q2 = f"Sum([EcotypeFrac]*(([CropFrac]*[CropEcoLDF{j}])\
  !        +([TreeFrac]*[TreeEcoLDF{j}])\
  !        +([HerbFrac]*[HerbEcoLDF{j}])\
  !        +([ShrubFrac]*[ShrubEcoLDF{j}]))) AS LDF{j}, "
  !        query_str += q2
  !    query_str = query_str.rstrip(', ')
  !    query_str += " FROM ((((GridGrowthForm INNER JOIN GridEcotype ON GridGrowthForm.gridID = GridEcotype.gridID)  \
  !    INNER JOIN IntermediateHerbEcoEF ON GridEcotype.EcotypeID = IntermediateHerbEcoEF.EcoTypeID)                  \
  !    INNER JOIN IntermediateShrubEcoEF ON IntermediateHerbEcoEF.EcoTypeID = IntermediateShrubEcoEF.EcoTypeID)      \
  !    INNER JOIN IntermediateTreeEcoEF ON IntermediateShrubEcoEF.EcoTypeID = IntermediateTreeEcoEF.EcoTypeID)       \
  !    INNER JOIN IntermediateCropEcoEF ON IntermediateTreeEcoEF.EcoTypeID = IntermediateCropEcoEF.EcoTypeID         \
  !    GROUP BY GridGrowthForm.gridID;"
  !
  !===============================================================================================================
  subroutine build_EFS_LDF(g,p,ecotypefile,cropfile,treefile,grassfile,shrubfile)
     implicit none
     type(grid_type) ,intent(in) :: g
     type(proj_type) ,intent(in) :: p
     character(len=200) :: laivfile
     real, allocatable :: DATA(:,:,:)
     integer :: date_time_dim_id,var_dim_id,lev_id,tstep_dim_id,var_id,ncid
     character(len=25),allocatable :: var_list(:),var_desc(:),var_unit(:)
     character(len=10) :: outfile
     integer :: nvars
     character(len=2)::kk
     real :: missing_value
     
     character(len=50) :: speciation_filename(num_growthforms)
     real, dimension(num_ef)  :: intermediate_ef_herb , intermediate_ef_shrub,  intermediate_ef_tree,  intermediate_ef_crop
     real, dimension(num_ldf) :: intermediate_ldf_herb, intermediate_ldf_shrub, intermediate_ldf_tree, intermediate_ldf_crop
     
     print*,"Building MEGAN_EFS & MEGAN_LDF ..."
     
      outfileEF='EFMAP.nc'
     outfileLDF='LDF.nc'

     !nvars=19  !nvars=4
     
     allocate(var_list(nvars))  
     allocate(var_unit(nvars))  
     allocate(var_desc(nvars))  
     allocate(LAIv(grid%nx,grid%ny,nvars))  
                                                                             
     file_list=(/ecotypefile,cropfile,treefile,grassfile,shrubfile/)
     !Levanto netcdf input files
     do k=1,nvars
         write(kk,'(I0.2)') k
         call check(nf90_open(trim(laivfile)//kk//".nc", nf90_write, ncid ))
              call check(   nf90_inq_varid(ncid,'Band1', var_id ))
              call check(   nf90_get_var(ncid, var_id , LAIv(:,:,k)  ))
              where (LAIv < 0.0 )
                      LAIv=0.0
              endwhere
         call check(nf90_close(ncid ))
     enddo


     !Bases de datos:
                     EF_file= "db/EFv210806.csv"             !Archivo con EF y LDF para cada VegId
     speciation_filename(1) = "db/SpeciationCrop210806.csv"  !"crop_speciation.csv"
     speciation_filename(2) = "db/SpeciationHerb210806.csv"  !"herb_speciation.csv"
     speciation_filename(3) = "db/SpeciationShrub210806.csv" !"shrub_speciation.csv"
     speciation_filename(4) = "db/SpeciationTree210725.csv"  !"tree_speciation.csv"
 

     !PSEUDO_CODE:=========================================================
     !para EF:    
     do i=ef_s, ef_e:   !para "i" en cada especie
         !<GF>_EcoEF(i) = sum( EF(i) * growthformFile%Specfrac )
         CropEcoEF_EcoEF(i)=sum( EF(i) * growthformFile%Specfrac )
         TreeEcoEF_EcoEF(i)=sum( EF(i) * growthformFile%Specfrac )
         HerbEcoEF_EcoEF(i)=sum( EF(i) * growthformFile%Specfrac )
        ShrubEcoEF_EcoEF(i)=sum( EF(i) * growthformFile%Specfrac )
      !(Join EF on VegID), Group by EcotypeId
      enddo

      !para LDF:                                                            
      do j=ldf_s,ldf_e  !para "j" en cada especie
         !<GF>_EcoLDF(j)=sum( LDF(j) * growthformFile%Specfrac )
         CropEcoEF_EcoLDF(j)=sum( LDF(j) * growthformFile%Specfrac )
         TreeEcoEF_EcoLDF(j)=sum( LDF(j) * growthformFile%Specfrac )
         HerbEcoEF_EcoLDF(j)=sum( LDF(j) * growthformFile%Specfrac )
        ShrubEcoEF_EcoLDF(j)=sum( LDF(j) * growthformFile%Specfrac )
      !(Join LDF on VegID), Group by EcotypeId
      enddo

     do i=ef_s, ef_e: !para "i" en cada especie
       EF(i) = Sum( [EcotypeFrac] *
                      (   [CropFrac ]*  CropEcoEF{i}  
                        + [TreeFrac ]*  TreeEcoEF{i} 
                        + [HerbFrac ]*  HerbEcoEF{i} 
                        + [ShrubFrac]* ShrubEcoEF{i} 
                       )
       )
     enddo
     do j=ldf_s, ldf_e: !para "j" en cada especie
       EF(i) = Sum( [EcotypeFrac] *
                      (   [CropFrac ]*  CropEcoEF{i}  
                        + [TreeFrac ]*  TreeEcoEF{i} 
                        + [HerbFrac ]*  HerbEcoEF{i} 
                        + [ShrubFrac]* ShrubEcoEF{i} 
                       )
       )
     enddo
     !END PSEUDO_CODE====================================================

     !Idea de ChatGPT:
     ! Calculate intermediate EF values for each growth form and ecotype
     !DO k = 1, num_growthforms
     !  OPEN(unit=k, file=speciation_filename(k))
     !  DO i = 1, num_ecotypes
     !    READ(k, *) ecotype_id, veg_id
     !    READ(k, *) intermediate_ef_crop(i)
     !    READ(k, *) intermediate_ef_herb(i)
     !    READ(k, *) intermediate_ef_shrub(i)
     !    READ(k, *) intermediate_ef_tree(i)
     !    READ(k, *) intermediate_ldf_crop(i)
     !    READ(k, *) intermediate_ldf_herb(i)
     !    READ(k, *) intermediate_ldf_shrub(i)
     !    READ(k, *) intermediate_ldf_tree(i)
     !  END DO
     !  CLOSE(k)
     !  
     !! Calculate final EF and LDF values for this growth form
     !final_ef = 0.0
     !final_ldf = 0.0
     !DO i = 1, num_ecotypes
     !  final_ef = final_ef + ecotype_frac(i) * &
     !      (growthform_frac(1) * intermediate_ef_crop(i) + &
     !      growthform_frac(2) * intermediate_ef_herb(i) + &
     !      growthform_frac(3) * intermediate_ef_shrub(i) + &
     !      growthform_frac(4) * intermediate_ef_tree(i))
     !  final_ldf = final_ldf + ecotype_frac(i) * &
     !      (growthform_frac(1) * intermediate_ldf_crop(i) + &
     !      growthform_frac(2) * intermediate_ldf_herb(i) + &
     !      growthform_frac(3) * intermediate_ldf_shrub(i) + &
     !      growthform_frac(4) * intermediate_ldf_tree(i))
     !END DO
     !
     !! Save final EF and LDF values for this growth form
     !DO i = 1, num_ef
     !  grid_ef(i) = final_ef(i)
     !END DO


    var_list=(/"EF_ISOP","EF_MBO","EF_MT_PINE","EF_MT_ACYC","EF_MT_CAMP","EF_MT_SABI","EF_MT_AROM","EF_NO","EF_SQT_HR","EF_SQT_LR","EF_MEOH","EF_ACTO","EF_ETOH","EF_ACID","EF_LVOC","EF_OXPROD","EF_STRESS","EF_OTHER","EF_CO"/)
    var_list=(/"EF_ISOP","EF_MBO","EF_MT_PINE","EF_MT_ACYC","EF_MT_CAMP","EF_MT_SABI","EF_MT_AROM","EF_NO","EF_SQT_HR","EF_SQT_LR","EF_MEOH","EF_ACTO","EF_ETOH","EF_ACID","EF_LVOC","EF_OXPROD","EF_STRESS","EF_OTHER","EF_CO"/)

     var_unit=(/"nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension","nondmension" /)
  
     ! Create the NetCDF file
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
  
      call check(nf90_close( ncid ))
      !Cierro NetCDF outFile
  
  end subroutine



























end program


!ESTRUCTURA GENERAL de MEGAN_PREP_Code_Jan_2022:
!(1) Leen Namelist
!(2) call wrf_file         (leen wrfout para sacar algunos parametros de grilla y proyeccion)
!(3) call  megan2_bioemiss (creo que aca es donde interpola)
!(4) call write_*
!       prepmegan4cmaq_cantype.f90:      write(99,'(a)')"CID,ICELL,JCELl,NEEDL,TROPI,BROAD,SHRUB,HERB,CROP"
!       prepmegan4cmaq_lai.f90:          write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,LAI01,LAI02,LAI03,LAI04,LAI05,LAI06,LAI07,LAI08,LAI09,LAI10, &
!       prepmegan4cmaq_ecotype.f90:      write(99,'(a)')"gridID,EcotypeID,EcotypeFrac"
!       prepmegan4cmaq_grwform.f90:      write(99,'(a)')"gridID,TreeFrac,CropFrac,ShrubFrac,HerbFrac"
!       prepmegan4cmaq_ef.f90:           write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,DSRAD,DTEMP,ISOP,MYRC,SABI,LIMO,A_3CAR,OCIM,BPIN,APIN,OMTP,FARN,BCAR,OSQT,MBO,MEOH,ACTO,CO,NO,BIDER,STRESS,OTHER"
!       prepmegan4cmaq_fert.f90:         write(99,888)"CELL_ID,X,Y,LAT,LONG",fertday
!       prepmegan4cmaq_landtype.f90:     write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,LANDTYP"
!       prepmegan4cmaq_nitrogen.f90:     write(99,'(a)')"CELL_ID,X,Y,LAT,LONG,NITROGEN01,NITROGEN02,NITROGEN03,NITROGEN04,NITROGEN05,NITROGEN06,NITROGEN07,NITROGEN08,NITROGEN09,NITROGEN10, &
!       prepmegan4cmaq_arid.f90:         write(99,'(a)')"CID, ICELL, JCELL,  ARID"
!       prepmegan4cmaq_non_arid.f90:     write(99,'(a)')"CID, ICELL, JCELL,  NON_ARID"
!       prepmegan4cmaq_pft.f90:          write(99,'(a)')"CID, ICELL, JCELL,  NT_EG_TEMP,  NT_DC_BORL,  NT_EG_BORL, &
!       prepmegan4cmaq_w126.f90:         write(99,'(a)')"CID, ICELL, JCELL,  W126(ppm-hours)"
!(5) clean-up (deallocate)
