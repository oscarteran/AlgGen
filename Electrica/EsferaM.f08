!--------------------------------------------------------------------------------------------------!
! UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
! FACULTAD DE INGENIERÍA
! DIVISIÓN DE INGENIERÍA EN CIENCIAS DE LA TIERRA
! DEPARTAMENTO DE GEOFÍSICA
!
! Asignatura: Inversión de Datos Geofísicos
!   Profesor: Dr. Mauricio Nava Flores
!--------------------------------------------------------------------------------------------------!
! Función EsferaM (Nava-Flores, 2019).
!  - Basada en la subrutina dipole (Blakely, 1996).
!
! Calcula la intensidad de campo total causada por una esfera magnetizada uniformemente en un punto 
! de observación.
!
! La esfera está referenciada en un sistema rectangular derecho (X dirigido hacia el Norte y Z 
! hacia abajo).
!
! Parámetros de entrada:
!   - Coordenadas del punto de observación: x0,y0,z0
!   - Coordenadas del centro de la esfera: xc,yc,zc
!   - Radio de la esfera: a
!   - Ángulo de inclinación de la magnetización de la esfera: mi
!   - Ángulo de declinación de la magnetización de la esfera: md
!   - Intensidad de magnetización de la esfera: mg
!
! Unidades consideradas:
!        - Longitud: Son irrelevantes pero deben ser consistentes
!        - Ángulos: [Grados]
!        - Intensidad de magnetización: [A/m]
!        - Inducción magnética: [nT]
!
! Últimas modificaciones: Nava-Flores, 23/Ene/2021
!--------------------------------------------------------------------------------------------------!
function EsferaM(x0, y0, z0, xc, yc, zc, a, mi, md, mg)

   use NumKind

   implicit none
   real(dp), intent(in):: x0, y0, z0, xc, yc, zc, a, mi, md, mg
   real(dp):: rx, ry, rz, r2, r, r5, dot, bx, by, bz, xmi, xmd, mx, my, mz, moment
   real(dp), parameter:: pi=4._dp*atan(1._dp), d2rad=pi/180._dp, t2nt=1.e9_dp
   real(dp):: EsferaM

   EsferaM=0._dp
   
   xmi=mi*d2rad
   xmd=(90._dp-md)*d2rad 

   mx=cos(xmi)*cos(xmd)
   my=cos(xmi)*sin(xmd)
   mz=sin(xmi)
   
   rx=x0-xc
   ry=y0-yc
   rz=z0-zc
   r2=rx**2+ry**2+rz**2
   r=sqrt(r2)

   r5=r**5
   dot=rx*mx+ry*my+rz*mz
   
   moment=4._dp*pi*(a**3)*mg/3._dp
   
   bx=(1.e-7_dp*moment*(3._dp*dot*rx-r2*mx)/r5)*t2nt
   by=(1.e-7_dp*moment*(3._dp*dot*ry-r2*my)/r5)*t2nt
   bz=(1.e-7_dp*moment*(3._dp*dot*rz-r2*mz)/r5)*t2nt
   
   EsferaM=bx*mx+by*my+bz*mz
   
end function EsferaM
!--------------------------------------------------------------------------------------------------!
