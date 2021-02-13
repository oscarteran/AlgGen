!--------------------------------------------------------------------------------------------------!
! UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
! FACULTAD DE INGENIERÍA
! DIVISIÓN DE INGENIERÍA EN CIENCIAS DE LA TIERRA
! DEPARTAMENTO DE GEOFÍSICA
!
! Asignatura: Inversión de Datos Geofísicos
!   Profesor: Dr. Mauricio Nava Flores
!--------------------------------------------------------------------------------------------------!
! Función Función CilindroG.f08 (Nava-Flores, 2018).
! 
! Útil para calcular la componente vertical de la aceleración gravitacional, debida a un cilindro 
! horizontal que se extiende infinitamente en dirección perpendicular al plano XZ, con densidad 
! constante, en un punto de observación.
!
! El cilindro está referenciado en un sistema rectangular derecho (Z hacia abajo).
!
! Parámetros de entrada:
!  - Puntos de observación: x0 (sólo coordenada x0; se asume y0=0 y z0=0)
!  - Parámetros que definen el cilindro horizontal:
!          r: Radio del cilindro
!         xc: Posición horizontal del eje del cilindro (c/r al vector de puntos de observación)
!         zc: Profundidad del eje del cilindro
!        rho: Contraste de densidad del cilindro con respecto al medio circundante
!
! Unidades consideradas:
!        - Longitud: [km]
!        - Densidad: [kg/m**3]
!        - Anomalía gravimétrica: [mGal]
!
! Últimas modificaciones: Nava-Flores, 23/Ene/2021
!--------------------------------------------------------------------------------------------------!
function CilindroG(x0, r, xc, zc, rho)

   use NumKind

   implicit none
   real(dp), intent(in):: x0, r, xc, zc, rho
   real(dp):: CilindroG
   real(dp), parameter:: pi=4._dp*atan(1._dp), gamma=6.674e-11_dp, si2mg=1.e5_dp, km2m=1.e3_dp

   CilindroG=2._dp*pi*gamma*r**2*rho*zc/((x0-xc)**2+zc**2)

   CilindroG=CilindroG*si2mg*km2m

end function CilindroG
!--------------------------------------------------------------------------------------------------!