module nc_handler_mod

        use netcdf
        use proj_mod
contains

 subroutine check(status)
   integer, intent(in) :: status
   if (status /= nf90_noerr) then
     write(*,*) nf90_strerror(status)
     stop 'netcdf error'
   end if
 end subroutine check

 subroutine createNetCDF(outFile,p,g,var_list,var_unit,var_desc)
    implicit none
    type(grid_type) , intent(in) :: g
    type(proj_type) , intent(in) :: p
    character(len=10), intent(in) :: outFile
    character(len=16), allocatable:: var_list(:),var_unit(:)
    character(len=25), allocatable:: var_desc(:)
    integer :: ncid,tstep_dim_id,date_time_dim_id,col_dim_id,row_dim_id,lay_dim_id,var_dim_id,var_id
    integer :: k
    integer :: nvars
    character(800) :: var_list_string
 
    nvars=size(var_list)
    write(var_list_string,*) var_list                !este es un global attr importante.
    
    call check(nf90_create(outFile, NF90_CLOBBER, ncid))
        ! Defino dimensiones
        call check(nf90_def_dim(ncid, "TSTEP"    , 1      , tstep_dim_id     ))
        call check(nf90_def_dim(ncid, "DATE-TIME", 2      , date_time_dim_id ))
        call check(nf90_def_dim(ncid, "COL"      , g%nx   , col_dim_id       ))
        call check(nf90_def_dim(ncid, "ROW"      , g%ny   , row_dim_id       ))
        call check(nf90_def_dim(ncid, "LAY"      , 1      , lay_dim_id       ))
        call check(nf90_def_dim(ncid, "VAR"      , nvars  , var_dim_id       ))
        !Defino variables
        call check(nf90_def_var(ncid,"TFLAG",NF90_INT      , [date_time_dim_id,var_dim_id,tstep_dim_id], var_id))
        call check(nf90_put_att(ncid, var_id, "units"      , "<YYYYDDD,HHMMSS>" ))
        call check(nf90_put_att(ncid, var_id, "long_name"  , "TFLAG           " ))
        call check(nf90_put_att(ncid, var_id, "var_desc"   , "Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS                                "))
        do k=1, nvars
          call check(nf90_def_var(ncid, trim(var_list(k)) , NF90_FLOAT, [col_dim_id,row_dim_id,lay_dim_id,tstep_dim_id], var_id)) 
          call check(nf90_put_att(ncid, var_id,"long_name",      var_list(k)  ))
          call check(nf90_put_att(ncid, var_id,"units"    , trim(var_unit(k)) ))
          call check(nf90_put_att(ncid, var_id,"var_desc" , trim(var_desc(k)) ))
        end do
        ! Defino attributos
        call check(nf90_put_att(ncid, nf90_global,"IOAPI_VERSION", "ioapi-3.2: \$Id: init3" ))
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
        call check(nf90_put_att(ncid, nf90_global,"NVARS"  , nvars                ))
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


end module
