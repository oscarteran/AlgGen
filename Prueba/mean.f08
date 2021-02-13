!--------------------------------------------------------------------------------------------------!
! UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
! FACULTAD DE INGENIERÍA
! DIVISIÓN DE INGENIERÍA EN CIENCIAS DE LA TIERRA
! DEPARTAMENTO DE GEOFÍSICA
!
! Asignatura: Inversión de Datos Geofísicos
!   Profesor: Dr. Mauricio Nava Flores
!--------------------------------------------------------------------------------------------------!
! Función para calcular la media aritmética de un conjunto de datos almacenado en un arreglo 1D.
!
! Últimas modificaciones: Nava-Flores, 20/Ene/2021
!--------------------------------------------------------------------------------------------------!
function mean(x)

   use NumKind

   implicit none
   real(dp), allocatable, intent(in):: x(:)
   real(dp):: mean

   mean=sum(x)/size(x,1)

end function mean
!--------------------------------------------------------------------------------------------------!