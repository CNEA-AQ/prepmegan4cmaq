!---------------------------------------------------------------
! author: Ramiro A. Espada. April 2023.
! Based on Prep_code &  MEGEFP32 (UCI-BAI-MEGAN)
!---------------------------------------------------------------
program prepmegan4cmaq

  use netcdf
  use utils_mod         !utils
  use PROJ_mod          !subroutines for coordinate transformations
  use nc_handler_mod    !functions to deal with netcdf
  !use INTERP_mod        !subroutines for interpolation/regridding

  implicit none

  type(proj_type) :: proj
  type(grid_type) :: grid

  integer :: status,iostat
  integer :: i,j,k

  character(len=17) :: start_date,end_date
  character(200)    :: griddesc_file,gridname
  character(200)    :: crop_frac_file,grass_frac_file,shrub_frac_file,tree_frac_file,nl_tree_frac_file,bl_tree_frac_file,tp_tree_frac_file,ecotype_file,laiv_files,GtEcoEF_file,arid_file, narid_file, lt_file,nitro_files

  namelist/control/start_date,end_date,griddesc_file,gridname,crop_frac_file,grass_frac_file,shrub_frac_file,tree_frac_file,nl_tree_frac_file, bl_tree_frac_file,tp_tree_frac_file,ecotype_file,laiv_files,GtEcoEF_file,arid_file, narid_file, lt_file,nitro_files

  !Leo namelist:
  read(*,nml=control, iostat=iostat)
  if( iostat /= 0 ) then
    write(*,*) 'prepmegan4cmaq: failed to read namelist; error = ',iostat
    stop
  end if

  !Leo GRIDDESC:
  call read_GRIDDESC(griddesc_file,gridname, proj, grid)

  call set_additional_proj_params(proj)
  call set_additional_grid_params(proj, grid)
                                                                       
  !`MEGAN_CTS` (*Canopy Type Fractions*) 
  call build_CT3(grid,proj,crop_frac_file,tree_frac_file,grass_frac_file,shrub_frac_file,nl_tree_frac_file,tp_tree_frac_file)
  
  !`MEGAN_LAI` (Leaf Area Index).
  !call build_LAIv(grid,proj,laiv_files)

  !`MEGAN_EFS` (emission factors) & `MEGAN_LDF` (*Light Dependence Fractions*) 
  call build_EFS_LDF(grid,proj,GtEcoEF_file,ecotype_file,crop_frac_file,tree_frac_file,grass_frac_file,shrub_frac_file)
  
  !BDSNP:
  !call BDSNP_AFILE()   !int Arid (0/1) &  call BDSNP_NAFILE()  !int Non-Arid (0/1) &  call BDSNP_LFILE()   !int Land types (1:24)
  call BDSNP_LAND(grid,proj,arid_file,narid_file,lt_file)

  call BDSNP_NFILE(grid,proj,nitro_files)   !float nitrogeno01, nitrogeno02,...,nitrogeno12  monthly nitrogen deposition in ng of N /m2/s
  !call BDSNP_FFILE()                       !float fert01,fert02,...,fert  daily fertilizer aplication. unit: ng of N/m2 

print*, "========================================="
print*, " prepmegan4cmaq: Completed successfully"
print*, "========================================="

contains

 subroutine read_GRIDDESC(griddescFile,gridName, p, g)
  implicit none
  character(200),intent(in) :: griddescFile
  character(*) ,intent(in)  :: gridName
  type(proj_type), intent(inout) :: p
  type(grid_type), intent(inout) :: g
  character(20) :: row
  iostat=0
  open(unit=2,file=griddescFile,status='old',action='read',access='sequential')
  do while(iostat == 0)  !loop por cada fila
     read(2,*,iostat=iostat) row
     if ( trim(row) == trim(gridname)) then
       g%gName=row
       read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%nx,g%ny !projName xorig yorig xcell ycell nrows ncols
       rewind(2)
     endif
     if (trim(row) == trim(p%pName)) then
       read(2,*) p%typ,p%alp,p%bet,p%gam,p%xcent,p%ycent   !map_proj truelat1 truelat2 stand_lon ref_lon ref_lat
       iostat=1
     endif
  enddo
  close(2)
 end subroutine

 function interpolate(p,g,inp_file)   result (img)
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(*),intent(in) :: inp_file
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
 !----------------------------------
 !  MEGAN_CTS 
 !---------------------------------
 subroutine build_CT3(g,p,cropfile,treefile,grassfile,shrubfile,nltreefile,troptreefile)
    implicit none
    type(grid_type) ,intent(in) :: g
    type(proj_type) ,intent(in) :: p
    character(len=200) :: cropfile,treefile,grassfile,shrubfile,nltreefile,troptreefile
    real, allocatable :: CTS(:,:,:,:)!,CTS(:,:,:)
    character(len=10) :: outfile
    character(len=16),allocatable :: var_list(:),var_unit(:)
    character(len=25),allocatable :: var_desc(:)
    integer :: nvars
    !integer :: var_id,ncid
    integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,var_id
 
    print*,"Building MEGAN_CTS file ..."

    outfile='CT3.nc'
    nvars=6
 
    allocate(var_list(nvars))  
    allocate(var_unit(nvars))  
    allocate(var_desc(nvars))  
    allocate(CTS(g%nx,g%ny,1,nvars+1))  ! allocate(CTS(g%nx,g%ny,nvars))  
  
    CTS(:,:,1,1)=interpolate(p,g,inp_file="./input/GF3aTree.nc"            ) !"GF3aTree.nc")   !,varname)
    CTS(:,:,1,2)=interpolate(p,g,inp_file="./input/GF3aShrub.nc"           ) !   shrubfile)   !,varname)
    CTS(:,:,1,3)=interpolate(p,g,inp_file="./input/GF3aCrop.nc"            ) !    cropfile)   !,varname)
    CTS(:,:,1,4)=interpolate(p,g,inp_file="./input/GF3aGrass.nc"           ) !   grassfile)   !,varname)
    CTS(:,:,1,5)=interpolate(p,g,inp_file="./input/NTfrac_reorder_lat.nc"  ) !  nltreefile)   !,varname)
    CTS(:,:,1,6)=interpolate(p,g,inp_file="./input/tropfrac_reorder_lat.nc") !troptreefile)   !,varname)

    !needleleaf tree
    CTS(:,:,1,5)=CTS(:,:,1,1) * (1.0-CTS(:,:,1,6)) * CTS(:,:,1,5)    
    !tropical tree
    CTS(:,:,1,6)=CTS(:,:,1,1) * CTS(:,:,1,6) 
    !boradleaf tree
    CTS(:,:,1,7)=CTS(:,:,1,1) * (1.0-CTS(:,:,1,6)) * (1.0-CTS(:,:,1,5))
    where ( CTS < 0.0 )
         CTS=0.0
    end where
    var_list=(/ 'SHRUB      ','CROP       ','HERB       ','NEEDL      ','TROPI      ','BROAD      '/)
    var_desc=(/ "shrub fraction         ","crop fraction          ","grass fraction         ", "needle tree fraction   ","tropical tree fraction ","broadleaf tree fraction" /) 
    var_unit=spread("nondimension",1,nvars)
 
    ! Create the NetCDF file
    call createNetCDF(outFile,p,g,var_list,var_unit,var_desc)

    !Abro NetCDF outFile
    call check(nf90_open(outFile, nf90_write, ncid       ))
    do i=2,nvars+1
       call check(nf90_inq_varid(ncid,TRIM(var_list(i-1)),var_id))
       call check(nf90_put_var(ncid, var_id, CTS(:,:,1,i)  ))        
    end do

    !TFLAG:
    call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
    call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,nvars) ))
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
    character(len=16),allocatable :: var_list(:),var_unit(:)
    character(len=25),allocatable :: var_desc(:)
    character(len=10) :: outfile
    integer ::var_id,ncid
    integer :: nvars
    character(len=2):: kk
    
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
        print*,trim(laivfile)//kk//".nc"
        LAIv(:,:,k)=interpolate(p,g,inp_file="laiv2003"//kk//"_30sec.nc") !"GF3aTree.nc")   !,varname)
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
      character(len=16), allocatable :: var_list(:),var_unit(:) !19 EF + 4 LDF
      character(len=25), allocatable :: var_desc(:) !19 EF + 4 LDF
      integer :: nvars
      real  :: EF(23)

      print*,"Building MEGAN_EFS & MEGAN_LDF ..."
      
      nvars=23
       outfileEF='EFMAP.nc'
      outfileLDF='LDF.nc'
 
      allocate( ECOTYPE(g%nx, g%ny   ))   !ecotype         
      allocate(    GTYP(g%nx, g%ny,4 ))   !growthtype fracs

      ECOTYPE(:,:)=interpolate(p,g, trim(ecotypefile))

      GTYP_LIST=(/'crop ','tree ','herb ','shrub'/)
      GTYP(:,:,1)=interpolate(p,g, cropfile)
      GTYP(:,:,2)=interpolate(p,g, treefile)
      GTYP(:,:,3)=interpolate(p,g,grassfile)
      GTYP(:,:,4)=interpolate(p,g,shrubfile)
      
      where (GTYP < 0.0 )
              GTYP=0.0
      endwhere
      GTYP=GTYP/100  !las frac están en porcentaje.

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
   
     outFileEF='EF.nc'
     ! Create the NetCDF file
     call createNetCDF(outFileEF,p,g,var_list,var_unit,var_desc)
     !Abro NetCDF outFile
     call check(nf90_open(outFileEF, nf90_write, ncid       ))
      do k=1, 23       
        call check(nf90_inq_varid(ncid,var_list(k) ,var_id ))
        call check(nf90_put_var(ncid, var_id, OUTGRID(:,:,k)))
      enddo
      !TFLAG:
      call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
      call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,19) ))
     call check(nf90_close( ncid ))
     !Cierro NetCDF outFile================

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
     character(len=16),allocatable :: var_list(:),var_unit(:)
     character(len=25),allocatable :: var_desc(:)
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

     LANDGRID(:,:,1) = interpolate(p,g,  arid_file)
     LANDGRID(:,:,2) = interpolate(p,g, narid_file)
     
     LANDGRID(:,:,3) =0.0 
     do lt=1,24
             write(landtypeId, '(I0.2)') lt
             LANDGRID(:,:,3) = LANDGRID(:,:,3)+lt*interpolate(p,g, trim(lt_file)//landtypeId//".nc") !DO AN INTEGER VERSION OF THIS
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
    character(len=16),allocatable :: var_list(:),var_unit(:)
    character(len=25),allocatable :: var_desc(:)
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
        !call check(nf90_open(trim(nitro_files)//kk//".nc", nf90_write, ncid ))
        !     call check(   nf90_inq_varid(ncid,'Band1', var_id ))
        !     call check(   nf90_get_var(ncid, var_id , NITRO(:,:,k)  ))
        !call check(nf90_close(ncid ))
        NITRO(:,:,k)  = interpolate(p,g,trim(nitro_files)//kk//".nc")
    
    enddo
    where (NITRO < 0.0 )
            NITRO=0.0
    endwhere
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
    !TFLAG:
    call check(nf90_inq_varid(ncid, "TFLAG"    , var_id))
    call check(nf90_put_var(ncid, var_id, spread((/0000000,000000 /),2,nvars) ))
    call check(nf90_close( ncid ))

  end subroutine 
end program
