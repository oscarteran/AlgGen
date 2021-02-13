!--------------------------------------------------------------------------------------------------!
! UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
! FACULTAD DE INGENIERÍA
! DIVISIÓN DE INGENIERÍA EN CIENCIAS DE LA TIERRA
! DEPARTAMENTO DE GEOFÍSICA
!
! Asignatura: Inversión de Datos Geofísicos
!   Profesor: Dr. Mauricio Nava Flores
!--------------------------------------------------------------------------------------------------!
! Función para convertir números binarios (cadena de caracteres) a enteros.
!
! Últimas modificaciones: Nava-Flores, 20/Ene/2021
!--------------------------------------------------------------------------------------------------!
function bin2dec(cadenabin)

   use NumKind

   implicit none
   character(len=*), intent(in):: cadenabin
   integer(il):: bin2dec

   read(cadenabin,'(b200)') bin2dec
   
end function bin2dec
!--------------------------------------------------------------------------------------------------!