!---------------------------------------------------------------
! author: Ramiro A. Espada. April 2023.
! Based on Prep_code &  MEGEFP32 (UCI-BAI-MEGAN)
!---------------------------------------------------------------
program prepmegan4cmaq

  use netcdf
  use utils_mod         !utils
  use proj_mod          !coordinate transformation functions
  use nc_handler_mod    !functions to deal with netcdf files
  use interpolate_mod   !function for interpolation/regridding
  use readGRIDDESC_mod  !GRIDDESC reader

  implicit none

  type(proj_type) :: proj
  type(grid_type) :: grid

  integer :: status,iostat
  integer :: i,j,k

  !character(len=17) :: start_date,end_date
  character(200)    :: griddesc_file,gridname,ecotypes_file,growtype_file,laiv_file,climate_file,fert_file,landtype_file,nitro_file,GtEcoEF_file
  logical           :: run_CTS=.true.,run_LAI=.true.,run_EFS=.true.,run_BDSNP=.true.
  character(6)      :: out_fmt='netcdf'
  namelist/control/griddesc_file,gridname,ecotypes_file,growtype_file,laiv_file,climate_file,landtype_file,nitro_file,fert_file,GtEcoEF_file,run_CTS,run_LAI,run_EFS,run_BDSNP,out_fmt

  !Leo namelist:
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'prepmegan4cmaq: failed to read namelist; error = ',iostat
    stop
  end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file,gridname, proj, grid)
                                                                       
  !`MEGAN_CTS` (*Canopy Type Fractions*) 
  if (run_CTS) call build_CT3(grid,proj,growtype_file)
  
  !`MEGAN_LAI` (Leaf Area Index).
  if (run_LAI) call build_LAIv(grid,proj,laiv_file)

  !`MEGAN_EFS` (Emission Factors) & `MEGAN_LDF` (*Light Dependence Factirs*) 
  if (run_EFS ) call build_EFS_LDF(grid,proj,GtEcoEF_file,ecotypes_file,growtype_file)
  
  !BDSNP: !call BDSNP_AFILE() !int Arid (0/1) &  call BDSNP_NAFILE()  !int Non-Arid (0/1) &  call BDSNP_LFILE()  !int Land types (1:24)
  if (run_BDSNP) then
     call BDSNP_LAND(grid,proj,climate_file,landtype_file)
     call BDSNP_NFILE(grid,proj,nitro_file)                !float nitrogeno01, nitrogeno02,...,nitrogeno12  monthly nitrogen deposition in ng of N /m2/s
     call BDSNP_FFILE(grid,proj,fert_file)                 !float fert01,fert02,...,fert  daily fertilizer aplication. unit: ng of N/m2 
  endif

print*, "========================================="
print*, " prepmegan4cmaq: Completed successfully"
print*, "========================================="

contains

 !----------------------------------
 !  MEGAN_CTS 
 !---------------------------------
 subroutine build_CT3(g,p,gwt_file) 
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(len=200) :: gwt_file 
    real, allocatable :: CTS(:,:,:,:)
    character(len=16) :: outfile
    character(len=16),allocatable :: var_list(:),var_unit(:)
    character(len=25),allocatable :: var_desc(:)
    integer :: nvars
    integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,var_id
 
    print*,"Building MEGAN_CTS file ..."

    nvars=6
    
    allocate(CTS(g%nx,g%ny,1,nvars+1))  ! allocate(CTS(g%nx,g%ny,nvars))  

    CTS(:,:,1,7)=interpolate(p,g,inp_file=gwt_file, varname="tree"     , method="bilinear") 
    CTS(:,:,1,1)=interpolate(p,g,inp_file=gwt_file, varname="nl_tree"  , method="bilinear")
    CTS(:,:,1,2)=interpolate(p,g,inp_file=gwt_file, varname="trop_tree", method="bilinear")
    CTS(:,:,1,4)=interpolate(p,g,inp_file=gwt_file, varname="shrub"    , method="bilinear")
    CTS(:,:,1,5)=interpolate(p,g,inp_file=gwt_file, varname="grass"    , method="bilinear")
    CTS(:,:,1,6)=interpolate(p,g,inp_file=gwt_file, varname="crop"     , method="bilinear")
    !where ( CTS < 0 .or. CTS > 100 )
    !     CTS=0
    !end where
    
    !needleleaf tree
    CTS(:,:,1,1)=(    CTS(:,:,1,1)/100.0) * CTS(:,:,1,7) * (1.0-CTS(:,:,1,2)/100.0) 
    !boradleaf tree
    CTS(:,:,1,3)=CTS(:,:,1,7) * (1.0-CTS(:,:,1,2)/100.0) * (1.0-CTS(:,:,1,1)/100.0)
    !tropical tree
    CTS(:,:,1,2)=CTS(:,:,1,7) * (    CTS(:,:,1,2)/100.0)
    
    if ( trim(out_fmt) == 'csv' .or. trim(out_fmt) == 'CSV' ) then

        outfile='CT3.csv'
        OPEN(UNIT=1, FILE=outfile,STATUS="NEW",ACTION="WRITE")
        write(1,'(A)')"CID,ICELL,JCELl,NEEDL,TROPI,BROAD,SHRUB,HERB,CROP"
        k=0 !(cell_id)
        do j=1,g%ny
            do i=1,g%nx
               k=k+1
               write(1,'(3(I0,","), 6(F10.4,","))') k,i,j,CTS(i,j,1,1:6)
            enddo
        enddo
        CLOSE(UNIT=1)

    else if ( trim(out_fmt) == 'NetCDF' .or. trim(out_fmt) == 'netcdf' .or. trim(out_fmt) == 'NETCDF') then

         outfile='CT3.nc'
         !!"Improved" version:
         !allocate(var_list(nvars))  
         !allocate(var_unit(nvars))  
         !allocate(var_desc(nvars))  
         !var_list=(/ 'NEEDL      ','TROPI      ','BROAD      ','SHRUB      ','GRASS      ','CROP       '/)
         !var_desc=(/  "needle tree fraction   ","tropical tree fraction ","broadleaf tree fraction", "shrub fraction         ","grass fraction         ","crop fraction          "/) 
         !var_unit=spread("nondimension",1,nvars)
         !var_dtyp=spread("FLOAT",1,nvars)
 
         !! Create the NetCDF file
         !call createNetCDF(outFile,p,g,var_list,var_unit,var_desc,var_dtyp)

         !!Abro NetCDF outFile
         !call check(nf90_open(outFile, nf90_write, ncid       ))
         !do i=1,nvars
         !   call check(nf90_inq_varid(ncid,TRIM(var_list(i)),var_id))
         !   call check(nf90_put_var(ncid, var_id, CTS(:,:,1,i)  ))        
         !end do

         !!TFLAG:
         !call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
         !call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,nvars) ))
         !call check(nf90_close(ncid))
         !!Cierro NetCDF outFile
         !call check(nf90_close(ncid))
         !**********************************************************!
               !MEGAN OLDER VERSION:
               ! Create the NetCDF file
               call check(nf90_create(outFile, NF90_CLOBBER, ncid))
                   ! Defino dimensiones
                   call check(nf90_def_dim(ncid, "TSTEP"    , nvars  , tstep_dim_id     ))
                   call check(nf90_def_dim(ncid, "DATE-TIME", 2      , date_time_dim_id ))
                   call check(nf90_def_dim(ncid, "COL"      , g%nx   , col_dim_id       ))
                   call check(nf90_def_dim(ncid, "ROW"      , g%ny   , row_dim_id       ))
                   call check(nf90_def_dim(ncid, "LAY"      , 1      , lay_dim_id       ))
                   call check(nf90_def_dim(ncid, "VAR"      , 1      , var_dim_id       ))
                   !Defino variables
                   call check(nf90_def_var(ncid,"TFLAG",NF90_INT      , [date_time_dim_id,var_dim_id,tstep_dim_id], var_id))
                   call check(nf90_put_att(ncid, var_id, "units"      , "<YYYYDDD,HHMMSS>" ))
                   call check(nf90_put_att(ncid, var_id, "long_name"  , "TFLAG           " ))
                   call check(nf90_put_att(ncid, var_id, "var_desc"   , "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                "))
               
                   call check(nf90_def_var(ncid, 'CTS' , NF90_FLOAT, [col_dim_id,row_dim_id,lay_dim_id,tstep_dim_id], var_id))
                   call check(nf90_put_att(ncid, var_id,"long_name", "CTS" ))
                   call check(nf90_put_att(ncid, var_id,"units"    , "nondimension    " ))
                   call check(nf90_put_att(ncid, var_id,"var_desc" , "" ))
               
                   ! Defino attributos
                   call check(nf90_put_att(ncid, nf90_global,"IOAPI_VERSION", "ioapi-3.2: $Id: init3" ))
                   call check(nf90_put_att(ncid, nf90_global,"EXEC_ID", "????????????????"   ))
                   call check(nf90_put_att(ncid, nf90_global,"FTYPE"  , 1                    ))
                   call check(nf90_put_att(ncid, nf90_global,"SDATE"  , 0000000              ))!stat_date (int)
                   call check(nf90_put_att(ncid, nf90_global,"STIME"  , 000000               ))
                   call check(nf90_put_att(ncid, nf90_global,"WDATE"  , 2023001              ))
                   call check(nf90_put_att(ncid, nf90_global,"WTIME"  , 000000               ))
                   call check(nf90_put_att(ncid, nf90_global,"CDATE"  , 2023001              ))
                   call check(nf90_put_att(ncid, nf90_global,"CTIME"  , 000000               ))
                   call check(nf90_put_att(ncid, nf90_global,"TSTEP"  , 10000                ))
                   call check(nf90_put_att(ncid, nf90_global,"NTHIK"  , 1                    ))!no sé que es.
                   call check(nf90_put_att(ncid, nf90_global,"NCOLS"  , g%nx                 ))
                   call check(nf90_put_att(ncid, nf90_global,"NROWS"  , g%ny                 ))
                   call check(nf90_put_att(ncid, nf90_global,"NLAYS"  , 1                    ))!grid%nz
                   call check(nf90_put_att(ncid, nf90_global,"NVARS"  , 1                    ))!son 6, pero apiladas en la dimension temporal
                   call check(nf90_put_att(ncid, nf90_global,"GDTYP"  , p%typ                ))
                   call check(nf90_put_att(ncid, nf90_global,"P_ALP"  , p%alp                ))
                   call check(nf90_put_att(ncid, nf90_global,"P_BET"  , p%bet                ))
                   call check(nf90_put_att(ncid, nf90_global,"P_GAM"  , p%gam                ))
                   call check(nf90_put_att(ncid, nf90_global,"XCENT"  , p%xcent              ))
                   call check(nf90_put_att(ncid, nf90_global,"YCENT"  , p%ycent              ))
                   call check(nf90_put_att(ncid, nf90_global,"XORIG"  , g%xmin               ))
                   call check(nf90_put_att(ncid, nf90_global,"YORIG"  , g%ymin               ))
                   call check(nf90_put_att(ncid, nf90_global,"XCELL"  , g%dx                 ))
                   call check(nf90_put_att(ncid, nf90_global,"YCELL"  , g%dy                 ))
                   call check(nf90_put_att(ncid, nf90_global,"VGTYP"  , -9999                ))!no sé que es.
                   call check(nf90_put_att(ncid, nf90_global,"VGTOP"  , 0.                   ))!no sé que es.
                   call check(nf90_put_att(ncid, nf90_global,"VGLVLS" , [0., 0.]             ))!no sé que es.
                   call check(nf90_put_att(ncid, nf90_global,"GDNAM"  , g%gName              ))
                   call check(nf90_put_att(ncid, nf90_global,"UPNAM"  , "prepMegan4cmaq.exe" ))!no sé que es.
                   !call check(nf90_put_att_any(ncid, nf90_global,"VAR-LIST",nf90_char, 16, "CTS"))
                   call check(nf90_put_att(ncid, nf90_global,"VAR-LIST","CTS"))
                   call check(nf90_put_att(ncid, nf90_global,"FILEDESC" , "MEGAN input file"   ))
                   call check(nf90_put_att(ncid, nf90_global,"HISTORY"  , ""                   ))
               call check(nf90_enddef(ncid))
               !End NetCDF define mode
               call check(nf90_open(outFile, nf90_write, ncid       ))
                 !CTS
                 call check(nf90_inq_varid(ncid,"CTS" ,var_id))
                 call check(nf90_put_var(ncid, var_id, CTS(:,:,:,1:nvars) ))
                 !!TFLAG:
                 call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
                 call check(nf90_put_var(ncid, var_id, spread(spread((/0000000,000000/),2,nvars),2,1) ))
                 call check(nf90_put_var(ncid, var_id, reshape([0, 0, 0, 10000, 0, 20000, 0, 30000, 0, 40000, 0, 50000], [2,1,6])))
               call check(nf90_close(ncid))
         !**********************************************************!
    endif
    deallocate(CTS)            !Libero memoria
 
 end subroutine build_CT3
 
 !----------------------------------
 !  MEGAN_LAI 
 !---------------------------------
 subroutine build_LAIv(g,p,laiv_file)
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(len=200) :: laiv_file
    real, allocatable :: LAIv(:,:,:)
    character(len=16),allocatable :: var_list(:),var_unit(:),var_dtyp(:)
    character(len=25),allocatable :: var_desc(:)
    character(len=16) :: outfile
    integer ::var_id,ncid
    integer :: nvars
    character(len=10):: temporal_avg  !flag on LAIv global file that specifies temporal average of data
    character(len=2):: kk
    real :: lat,lon
    
    print*,"Building MEGAN_LAI file ..."
    
    !
    !Idea for implementing 8-day LAIv: just read global_var (NVARS) and depending of this read & write the output.
    call check(nf90_open(laiv_file,nf90_nowrite, ncid) )
    call check(nf90_get_att(ncid, NF90_GLOBAL, "temporal_average", temporal_avg) )
    call check(nf90_close(ncid))
    if ( trim(temporal_avg) == "monthly") then
       nvars=12
    else if ( trim(temporal_avg) == "8-day") then
       nvars=46
    else
       print*,"Couldn't get temporal_average (8-day or monthly) global attribute from LAI file",temporal_avg;stop
    endif
 
    allocate(var_list(nvars))  
    allocate(var_unit(nvars))  
    allocate(var_desc(nvars))  
    allocate(var_dtyp(nvars))  
    allocate(LAIv(g%nx,g%ny,nvars))  
 
    !Levanto netcdf input files
    do k=1,nvars
        write(var_list(k),'(A,I0.2)') "LAI",k
        write(var_desc(k),'(A,I0.2)') "LAI",k
        var_dtyp(k)="FLOAT"
        write(kk,'(I0.2)') k
        LAIv(:,:,k)=interpolate(p,g,inp_file=laiv_file,varname="laiv"//kk, method="bilinear") 
    enddo
    where (LAIv < 0.0 )
            LAIv=0.0
    endwhere
    var_unit=spread("nondimension",1,nvars)

    if ( trim(out_fmt) == 'csv' .or. trim(out_fmt) == 'CSV' ) then
        outfile='LAI3.csv'
        OPEN(UNIT=1, FILE=outfile,STATUS="NEW",ACTION="WRITE")
        write(1,'(A)') "CELL_ID,X,Y,LAT,LONG,LAI01,LAI02,..."
        k=0 !(cell_id)
        do j=1,g%ny
            do i=1,g%nx
                k=k+1
                call xy2ll(p, g%xmin+g%dx*i, g%ymin+g%dy*j,lon,lat)
                write(1,'(3(I10,","), 48(F10.4,","))') k,i,j,lat,lon,LAIv(i,j,:)/1000.0
            enddo
        enddo
        print*,'finished writing LAI to .csv file for CMAQ'
        CLOSE(UNIT=1)
                                                                                                            
    else if ( trim(out_fmt) == 'NetCDF' .or. trim(out_fmt) == 'netcdf' .or. trim(out_fmt) == 'NETCDF') then
        outfile='LAI3.nc'

        !Creo NetCDF file
        call createNetCDF(outFile,p,g,var_list,var_unit,var_desc,var_dtyp)
        
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
    endif
 end subroutine
 
  !----------------------------------
  !  MEGAN_EF & MEGAN_LDF
  !---------------------------------
   subroutine build_EFS_LDF(g,p,GtEcoEF_file,ecotypefile,gwt_file) !cropfile,treefile,grassfile,shrubfile)
      implicit none
      type(grid_type) ,intent(in) :: g
      type(proj_type) ,intent(in) :: p
      character(len=50),intent(in) :: ecotypefile,GtEcoEF_file, gwt_file !cropfile,treefile,grassfile,shrubfile
      character(len=16) :: outfileEF,outfileLDF,outfileEcotype, outfileGTfrac
      integer i,j,k
      integer :: var_id,ncid
      
      real,    allocatable :: GTYP(:,:,:)       !growth type fraction
      integer, allocatable :: ECOTYPE(:,:)      !este puede ser int tamb
      integer :: EcoID
      real, allocatable :: OUTGRID(:,:,:)

      character(len=5) :: GTYP_LIST(4) ! crop, tree, grass, shrub
      character(len=5) :: GtID
      character(len=16), allocatable :: var_list(:),var_unit(:),var_dtyp(:) !19 EF + 4 LDF
      character(len=25), allocatable :: var_desc(:) !19 EF + 4 LDF
      integer :: nvars
      real  :: EF(23)
      real ::lat,lon

      print*,"Building MEGAN_EFS & MEGAN_LDF ..."
      
      nvars=23
      
      allocate( ECOTYPE(g%nx, g%ny   ))   !ecotype         
      allocate(    GTYP(g%nx, g%ny,4 ))   !growthtype fracs

      ECOTYPE(:,:)=FLOOR(interpolate(p,g,ecotypefile,varname="ecotype", method="mode")) 

      GTYP_LIST=(/'crop ','tree ','grass','shrub'/)
      GTYP(:,:,1)=interpolate(p,g,gwt_file,varname="crop" , method="bilinear")
      GTYP(:,:,2)=interpolate(p,g,gwt_file,varname="tree" , method="bilinear")
      GTYP(:,:,3)=interpolate(p,g,gwt_file,varname="grass", method="bilinear")
      GTYP(:,:,4)=interpolate(p,g,gwt_file,varname="shrub", method="bilinear")
      
      where ( GTYP < 0.0 )
              GTYP=0.0
      endwhere
      GTYP=GTYP/100  ! % -> (0-1) . xq las frac están en porcentaje.

     !Arranco con una tabla: GtEcoEF_file:
     !growtypeId   ecotypeID   var1,   var2, ... ,  var19,  var20, ... ,  var23
     !crop         1           EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
     !crop         2           EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
     !crop         3           EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
     !....         ...         EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
     !tree         1700        EFi01, EFi02, ... ,  EFi19, LDFi01, ... , LDiF04
 
     allocate(OUTGRID(g%nx, g%ny, 23))   !outgrids EF1,EF2,...,LDF1,LDF2,..
     OUTGRID=0.0
     j=0
     open(unit=1,file=trim(GtEcoEF_file),status='unknown',action='read')
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

     allocate(var_list(23),var_desc(23),var_unit(23),var_dtyp(23))
     var_list=(/"EF_ISOP   ","EF_MBO    ","EF_MT_PINE","EF_MT_ACYC","EF_MT_CAMP","EF_MT_SABI","EF_MT_AROM","EF_NO     ","EF_SQT_HR ","EF_SQT_LR ", "EF_MEOH   ","EF_ACTO   ","EF_ETOH   ","EF_ACID   ","EF_LVOC   ","EF_OXPROD ","EF_STRESS ","EF_OTHER  ","EF_CO     ", "LDF03     ", "LDF04     ", "LDF05     ", "LDF06     " /)
     var_desc=var_list
     var_unit=spread("nanomol/m^2/s",1,nvars)
     var_dtyp=spread("FLOAT",1,nvars)
   
    if ( trim(out_fmt) == 'csv' .or. trim(out_fmt) == 'CSV' ) then
        outFileEF='EF.csv'
        outFileLDF='LDF.csv'
        outFileEcotype='ecotype.csv'
        outFileGTfrac='growth_form.csv'
        open(unit=1,file=outFileEF,status='unknown')
        open(unit=2,file=outFileLDF,status='unknown')
        open(unit=3,file=outFileEcotype,status='unknown')
        open(unit=4,file=outFileGTfrac,status='unknown')
        write(1,'(a)')"CELL_ID,X,Y,LAT,LONG,DSRAD,DTEMP,ISOP,MYRC,SABI,LIMO,A_3CAR,OCIM,BPIN,APIN,OMTP,FARN,BCAR,OSQT,MBO,MEOH,ACTO,CO,NO,BIDER,STRESS,OTHER"
        write(2,'(a)')"CELL_ID,X,Y,LAT,LONG,LDF01,LDF02,LDF03,LDF04"
        write(3,'(a)')"gridID,EcotypeID,EcotypeFrac"
        write(4,'(a)')"gridID,TreeFrac,CropFrac,ShrubFrac,HerbFrac"
        k=0 !(cell_id)
        do j=1,g%ny
            do i=1,g%nx
                k=k+1
                call xy2ll(p, g%xmin+g%dx*i, g%ymin+g%dy*j,lon,lat)
                write(1, '(3(I5,","), 2(F10.4,","),2(I0,","),19(F10.4,","),1(F10.4))') k,i,j,lat,lon,250,300, OUTGRID(i,j,1:19) !&
                   !MSEBIO_ISOP(ilon,ilat),MSEBIO_MYRC(ilon,ilat),MSEBIO_SABI(ilon,ilat), &
                   !MSEBIO_LIMO(ilon,ilat),MSEBIO_A_3CAR(ilon,ilat),MSEBIO_OCIM(ilon,ilat), &
                   !MSEBIO_BPIN(ilon,ilat),MSEBIO_APIN(ilon,ilat),MSEBIO_OMTP(ilon,ilat), &
                   !MSEBIO_FARN(ilon,ilat),MSEBIO_BCAR(ilon,ilat),MSEBIO_OSQT(ilon,ilat), &
                   !MSEBIO_MBO(ilon,ilat),MSEBIO_MEOH(ilon,ilat),MSEBIO_ACTO(ilon,ilat), &
                   !MSEBIO_CO(ilon,ilat),MSEBIO_NO(ilon,ilat),MSEBIO_BIDER(ilon,ilat), &
                   !MSEBIO_STRESS(ilon,ilat),MSEBIO_OTHER(ilon,ilat)
                write(2, '(3(I5,","), 3(F10.4,","),1(F10.4))') k,i,j,lat,lon, OUTGRID(i,j,20:23) 
                write(3, '(2(I0,","),F0.4)') k,ECOTYPE(i,j),1.0
                write(4, '(1(I0,","), 3(F0.4,","), F0.4)') k,GTYP(i,j,2),GTYP(i,j,1),GTYP(i,j,4),GTYP(i,j,3)
           end do
        end do
        print*,'jiang finished EF, LDF, ECOTYPE and GROWTYPE'
        close(1);close(2);close(3);close(4)

    else if ( trim(out_fmt) == 'NetCDF' .or. trim(out_fmt) == 'netcdf' .or. trim(out_fmt) == 'NETCDF') then
        outFileEF='EF.nc'
        ! Create the NetCDF file
        call createNetCDF(outFileEF,p,g,var_list(1:19),var_unit(1:19),var_desc(1:19),var_dtyp(1:19))
        !Abro NetCDF outFile
        call check(nf90_open(outFileEF, nf90_write, ncid       ))
         do k=1, 19       
           call check(nf90_inq_varid(ncid,var_list(k) ,var_id ))
           call check(nf90_put_var(ncid, var_id, OUTGRID(:,:,k)))
         enddo
         !TFLAG:
         call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
         call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,19) ))
        call check(nf90_close( ncid ))
        !Cierro NetCDF outFile================

        outFileLDF='LDF.nc'
        ! Create the NetCDF file
        call createNetCDF(outFileLDF,p,g,var_list(20:23),var_unit(20:23),var_desc(20:23),var_dtyp(20:23))
        !Abro NetCDF outFile
        call check(nf90_open(outFileLDF, nf90_write, ncid       ))
         do k=20, 23       
           call check(nf90_inq_varid(ncid,var_list(k) ,var_id ))
           call check(nf90_put_var(ncid, var_id, OUTGRID(:,:,k)))
         enddo
         !TFLAG:
         call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
         !call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,19) ))
         call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,4) ))
        call check(nf90_close( ncid ))
        !Cierro NetCDF outFile================
     endif
     deallocate(OUTGRID)
     deallocate(GTYP)
     deallocate(ECOTYPE)
 end subroutine

 !----------------------------------
 !  BDSNP_ARID, BDSNP_NONARID & BDSNP_LANDTYPE
 !---------------------------------
 subroutine BDSNP_LAND(g,p,climate_file,lt_file)
     implicit none
     type(grid_type) ,intent(in) :: g
     type(proj_type) ,intent(in) :: p
     character(len=200),intent(in) :: climate_file,lt_file
     real, allocatable :: LANDGRID(:,:,:)
     character(len=16),allocatable :: var_list(:),var_unit(:),var_dtyp(:)
     character(len=25),allocatable :: var_desc(:)
     integer :: var_id,ncid
     integer :: nvars
     integer :: lt
     character(len=2) ::landtypeId
     character(len=16) :: outfile, outfile_arid,outfile_nonarid,outfile_landtype
     real ::lat,lon

     print*,"Building BDSNP_ARID, BDSNP_NONARID & BDSNP_LANDTYPE ..."
     nvars=3
     allocate(LANDGRID(g%nx,g%ny,nvars))
     allocate(var_list(nvars))
     allocate(var_unit(nvars))
     allocate(var_desc(nvars))

     LANDGRID(:,:,3) = 1  
     LANDGRID(:,:,1) = interpolate(p,g,climate_file, varname="arid"    , method="mode")
     LANDGRID(:,:,2) = interpolate(p,g,climate_file, varname="non_arid", method="mode")
     LANDGRID(:,:,3) = interpolate(p,g,     lt_file, varname="landtype", method="mode")
     
     var_list=(/ 'ARID    ', 'NONARID ', 'LANDTYPE' /)
     var_unit=(/'1 or 0      ','1 or 0      ','nondimension' /) 
     var_desc=(/ 'ARID    ', 'NONARID ', 'LANDTYPE' /)
     var_dtyp=spread('INT',1,nvars) 
 
     if ( trim(out_fmt) == 'csv' .or. trim(out_fmt) == 'CSV' ) then
            outFile_arid='arid.csv'
        outFile_nonarid='non_arid.csv'
         outFile_landtype='landtype.csv'
        open(unit=1,file=outFile_arid    ,status='unknown')
        open(unit=2,file=outFile_nonarid ,status='unknown')
        open(unit=3,file=outFile_landtype,status='unknown')
        write(1,'(a)')"CID,ICELL,JCELL,ARID"
        write(2,'(a)')"CID,ICELL,JCELL,NON_ARID"
        write(3,'(a)')"CELL_ID,X,Y,LAT,LONG,LANDTYP"
        k=0 !(cell_id)
        do j=1,g%ny
            do i=1,g%nx
                k=k+1
                call xy2ll(p, g%xmin+g%dx*i, g%ymin+g%dy*j,lon,lat)
                write(1, '(3(I0,","), F10.4,",")') k,i,j,LANDGRID(i,j,1)
                write(2, '(3(I0,","), F10.4,",")') k,i,j,LANDGRID(i,j,2)
                write(3, '(3(I0,","), 2(F10.4,","),1(I0,","))') k,i,j,lat,lon,int(LANDGRID(i,j,3))
            enddo
        enddo
        close(1);close(2);close(3)

     else if ( trim(out_fmt) == 'NetCDF' .or. trim(out_fmt) == 'netcdf' .or. trim(out_fmt) == 'NETCDF') then
        outfile='LAND.nc'!(/"ARID.nc    ", "LANDTYPE.nc", "NONARID.nc "/) !
        ! Create the NetCDF file
        call createNetCDF(outFile,p,g,var_list,var_unit,var_desc,var_dtyp)

        !Abro NetCDF outFile
        call check(nf90_open(outFile, nf90_write, ncid       ))
        do k=1, nvars    
           call check(nf90_inq_varid(ncid,var_list(k) ,var_id ))
           call check(nf90_put_var(ncid, var_id, LANDGRID(:,:,k)))
        enddo
        !TFLAG:
        call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
        call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,nvars) ))
        call check(nf90_close( ncid ))
      endif
 end subroutine BDSNP_LAND

 subroutine BDSNP_NFILE(g,p,nitro_file)   
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(len=200),intent(in) :: nitro_file
    real, allocatable :: NITRO(:,:,:)
    character(len=16) :: outfile
    character(len=16),allocatable :: var_list(:),var_unit(:),var_dtyp(:)
    character(len=25),allocatable :: var_desc(:)
    integer :: var_id,ncid
    integer :: nvars
    integer :: k
    character(len=2) ::kk
    real ::lat,lon

    print*,"Building BDSNP_NFILE ..."
    nvars=12

    allocate(var_list(nvars))
    allocate(var_unit(nvars))
    allocate(var_desc(nvars))
    allocate(var_dtyp(nvars))
    allocate(NITRO(g%nx,g%ny,nvars))  
                                                                            
    !Levanto netcdf input files
    do k=1,nvars
        write(var_list(k),'(A,I0.2)') "NITROGEN",k
        write(var_desc(k),'(A,I0.2)') "NITROGEN",k
        var_dtyp(k)="FLOAT"
        write(kk,'(I0.2)') k
        NITRO(:,:,k)  = interpolate(p,g,nitro_file, varname="nitro"//kk, method="bilinear")
    enddo
    where (NITRO < 0.0 )
            NITRO=0.0
    endwhere
    NITRO=NITRO*1E+12! convert from kg/m2/s to ng/m2/s
    var_unit = spread( "ng/m2/s         ",1,nvars)
    var_desc=spread("monthly average total nitrogen deposition",1,nvars) 
    
    if ( trim(out_fmt) == 'csv' .or. trim(out_fmt) == 'CSV' ) then
      outFile='NITRO.csv'                              
      open(unit=1,file=outFile,status='unknown')
      write(1,'(a)')"CELL_ID,X,Y,LAT,LONG,NITROGEN01,NITROGEN02,NITROGEN03,NITROGEN04,NITROGEN05,NITROGEN06,NITROGEN07,NITROGEN08,NITROGEN09,NITROGEN10,NITROGEN11,NITROGEN12"
      k=0 !(cell_id)
      do j=1,g%ny
           do i=1,g%nx
               k=k+1
                call xy2ll(p, g%xmin+g%dx*i, g%ymin+g%dy*j,lon,lat)
               write(1, '(3(I0,","), 2(F10.4,","),12(E14.6,","))') k,i,j,lat,lon,NITRO(i,j,:)
          end do
      end do
      print*,'finished writing NITROGEN to .csv file for MEGAN'
      close(1)

    else if ( trim(out_fmt) == 'NetCDF' .or. trim(out_fmt) == 'netcdf' .or. trim(out_fmt) == 'NETCDF') then
      outfile='NDEP.nc'
      ! Create the NetCDF file
      call createNetCDF(outFile,p,g,var_list,var_unit,var_desc,var_dtyp)
                                                                 
      !Abro NetCDF outFile
      call check(nf90_open(outFile, nf90_write, ncid       ))
      do k=1, nvars    
        call check(nf90_inq_varid(ncid,var_list(k) ,var_id ))
        call check(nf90_put_var(ncid, var_id, NITRO(:,:,k)))
      enddo
      !TFLAG:
      call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
      call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,nvars) ))
      call check(nf90_close( ncid ))
    endif
  end subroutine 

  subroutine BDSNP_FFILE(g,p,fert_file)   
     implicit none
     type(grid_type) ,intent(in) :: g
     type(proj_type) ,intent(in) :: p
     character(len=200),intent(in) :: fert_file
     real, allocatable :: FERT(:,:,:)
     character(len=16) :: outfile
     character(len=16),allocatable :: var_list(:),var_unit(:),var_dtyp(:)
     character(len=25),allocatable :: var_desc(:)
     integer :: var_id,ncid
     integer :: nvars
     integer :: k
     character(len=3) ::kk
     real ::lat,lon
     print*,"Building BDSNP_FFILE ..."
     nvars=366
 
     allocate(var_list(nvars))
     allocate(var_unit(nvars))
     allocate(var_desc(nvars))
     allocate(var_dtyp(nvars))
     allocate(FERT(g%nx,g%ny,nvars))  
                                                                             
     !Levanto netcdf input files
     do k=1,nvars
         write(var_list(k),'(A,I0.3)') "FERT",k
         write(kk,'(I0.3)') k
         var_dtyp(k)="FLOAT"
         FERT(:,:,k)  = interpolate(p,g,fert_file, varname="fert"//kk, method="bilinear")
     enddo
     where (FERT< 0.0 )
            FERT=0.0
     endwhere
     !FERT=FERT*1E+6 !convert mg/m3 --> ng/m3 
     FERT=FERT*1E-6 !convert mg/m3 --> ng/m3 
     var_unit = spread( "ng/m2/s         ",1,nvars)
     var_desc=spread("monthly average total nitrogen deposition",1,nvars) 
     
     if ( trim(out_fmt) == 'csv' .or. trim(out_fmt) == 'CSV' ) then
       outFile='FERT.csv'                              
       open(unit=1,file=outFile,status='unknown')
       write(1,'(a)')"CELL_ID,X,Y,LAT,LONG,FERT001,FERT002,...,FERT366"
       k=0 !(cell_id)
       do j=1,g%ny
            do i=1,g%nx
                k=k+1
                call xy2ll(p, g%xmin+g%dx*i, g%ymin+g%dy*j,lon,lat)
                write (1, '(3(I0,","), 2(F10.4,","), 366(E10.4,","))') k,i,j,lat,lon,FERT(i,j,:)
           end do
       end do
       print*,'finished writing FERT to .csv file for MEGAN'
       close(1)
 
     else if ( trim(out_fmt) == 'NetCDF' .or. trim(out_fmt) == 'netcdf' .or. trim(out_fmt) == 'NETCDF') then

       outfile='FERT.nc'
       ! Create the NetCDF file
       call createNetCDF(outFile,p,g,var_list,var_unit,var_desc,var_dtyp)
                                                                  
       !Abro NetCDF outFile
       call check(nf90_open(outFile, nf90_write, ncid       ))
       do k=1, nvars    
         call check(nf90_inq_varid(ncid,var_list(k) ,var_id ))
         call check(nf90_put_var(ncid, var_id, FERT(:,:,k)))
       enddo
       !TFLAG:
       call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
       call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,nvars) ))
       call check(nf90_close( ncid ))
     endif
   end subroutine 

end program
