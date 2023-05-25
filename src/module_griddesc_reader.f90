module readGRIDDESC_mod
 use proj_mod
 implicit none

contains
 subroutine read_GRIDDESC(griddescFile,gridName, p, g)
  implicit none
  character(200),intent(in) :: griddescFile
  character(*) ,intent(in)  :: gridName
  type(proj_type), intent(inout) :: p
  type(grid_type), intent(inout) :: g
  character(20) :: row
  integer :: iostat
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
       read(2,*) p%typ,p%alp,p%bet,p%gam,p%xcent,p%ycent  
       iostat=1
     endif
  enddo
  close(2)

  call set_additional_proj_params(p)
  call set_additional_grid_params(p,g)
 end subroutine

 subroutine set_additional_proj_params(p)
    implicit none
    type(proj_type) ,intent(inout) :: p

    if ( p%typ == 2 ) then       !lambert conformal conic:        
       if ( ABS(p%alp - p%bet) > 0.1 ) then  !secant proj case
          p%p2=     LOG( COS(p%alp           *deg2rad )/ COS(p%bet           *deg2rad)   )
          p%p2=p%p2/LOG( TAN((45.0+0.5*p%bet)*deg2rad )/ TAN((45.0+0.5*p%alp)*deg2rad)   ) !n
        else                                 !tangent proj case
          p%p2=SIN(p%alp*deg2rad) !n
       endif
       p%p3=R_EARTH*(COS(p%alp*deg2rad)*TAN((45+0.5*p%alp)*deg2rad)**p%p2)*(1/p%p2)  !F
       p%p1=p%p3/(TAN((45 + 0.5*p%ycent)*deg2rad)**p%p2)                             !rho0 

    else if ( p%typ == 6 ) then  !polar secant stereographic
       print*, "Todavia no desarrollado soporte para proyeccion polar stereografica (ptype=",p%typ,")."; stop

    else if ( p%typ == 7 ) then  !equatorial mercator
       p%p1=COS(p%alp*deg2rad)   !k0

    else
        print*, "codigo de proyección invalido:",p%typ,"."; stop
    end if
 end subroutine

 subroutine set_additional_grid_params(p,g)
    implicit none
    type(proj_type) ,intent(inout) :: p
    type(grid_type) ,intent(inout) :: g
    real :: latmin,lonmin,latmax,lonmax
    !Obtener coordenadas del centro de la grilla, min y max:
    g%xc=0.0;g%yc=0.0;g%xmax=g%xmin+g%dx*g%nx; g%ymax=g%ymin+g%dy*g%ny

    !calculo minimos y maximos de latlon 
    !   (ojo! Dado que son transf no-lineales no corresponden necesariamente a los vertices)
    call xy2ll(p,g%xmin,g%ymin,g%lonmin,g%latmin)       !lower-left
    call xy2ll(p,g%xmax,g%ymax,g%lonmax,g%latmax)       !upper-right
 
    !latmin
    call xy2ll(p,g%xmin+g%dx*g%nx*0.5, g%ymin,lonmin,latmin)
    g%latmin=min(g%latmin,latmin)
    !latmax
    call xy2ll(p,g%xmax-g%dx*g%nx*0.5,g%ymax ,lonmax,latmax)
    g%latmax=max(g%latmax,latmax)

    !lonmin   
    call xy2ll(p,g%xmin,g%ymin+g%dy*g%ny*0.5,lonmin,latmin)
    g%lonmin=min(g%lonmin,lonmin)
    !         
    call xy2ll(p,g%xmin,g%ymax              ,lonmin,latmin)
    g%lonmin=min(g%lonmin,lonmin)

    !lonmax
    call xy2ll(p,g%xmax,g%ymax-g%dy*g%ny*0.5,lonmax,latmax)
    g%lonmax=max(g%lonmax,lonmax)
    !      
    call xy2ll(p,g%xmax,g%ymin              ,lonmax,latmax)
    g%lonmax=max(g%lonmax,lonmax)

 end subroutine

end module
