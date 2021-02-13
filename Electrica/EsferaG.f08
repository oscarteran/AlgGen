!--------------------------------------------------------------------------------------------------!
! UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
! FACULTAD DE INGENIERÍA
! DIVISIÓN DE INGENIERÍA EN CIENCIAS DE LA TIERRA
! DEPARTAMENTO DE GEOFÍSICA
!
! Asignatura: Inversión de Datos Geofísicos
!   Profesor: Dr. Mauricio Nava Flores
!--------------------------------------------------------------------------------------------------!
! Función EsferaG (Nava-Flores, 2017).
!
! Útil para calcular la componentes vertical de la aceleración gravitacional debida a una esfera con
! densidad constante, en un punto de observación.
!
! La esfera está referenciada en un sistema rectangular derecho (Z hacia abajo).
!
! Parámetros de entrada:
!  - Coordenadas del punto de observación: x0, y0, z0
!  -  Coordenadas del centro de la esfera: xc, yc, zc
!  -                   Radio de la esfera: a
!  -                Densidad de la esfera: Rho
!
! Unidades consideradas:
!        - Longitud: [km]
!        - Densidad: [kg/m**3]
!        - Vector de Atracción Gravitacional: [mGal]
!
! Últimas modificaciones: Nava-Flores, 23/Ene/2021
!--------------------------------------------------------------------------------------------------!
function EsferaG(x0, a, xc, zc,  Rho)

   use NumKind

   implicit none
   real(dp), intent(in):: x0, xc, zc, a, Rho
   real(dp):: k, R
   real(dp), parameter:: pi=4._dp*atan(1._dp), gamma=6.674e-11_dp, si2mg=1.e5_dp, km2m=1.e3_dp
   real(dp):: EsferaG

   k=(4._dp*pi*gamma*Rho*a**3)/3._dp

   R=sqrt((x0-xc)**2 + (zc)**2)

   EsferaG=(-k*(z0-zc)/(R**3))*si2mg*km2m
   
end function EsferaG
!--------------------------------------------------------------------------------------------------!
