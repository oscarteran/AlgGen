!--------------------------------------------------------------------------------------------------!
! UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
! FACULTAD DE INGENIERÍA
! DIVISIÓN DE INGENIERÍA EN CIENCIAS DE LA TIERRA
! DEPARTAMENTO DE GEOFÍSICA
!
! Asignatura: Inversión de Datos Geofísicos
!   Profesor: Dr. Mauricio Nava Flores
!--------------------------------------------------------------------------------------------------!
! Función para convertir números enteros a binarios (cadena de caracteres).
!
! Últimas modificaciones: Nava-Flores, 20/Ene/2021
!--------------------------------------------------------------------------------------------------!
function dec2bin(d,n)

   use NumKind

   implicit none
   integer(il), intent(in):: d
   integer(il), intent(in), optional:: n
   character(len=200):: bin
   integer(il):: long
   character(len=:), allocatable:: dec2bin

   write(bin,'(b0)') d
   long=len_trim(bin)

   if (present(n)) then
      do
         long=len_trim(bin)
         if (long >= n) exit
            bin='0'//trim(bin)
      end do
   end if

   allocate(character(long)::dec2bin)
   dec2bin=trim(bin)
   
end function dec2bin
!--------------------------------------------------------------------------------------------------!