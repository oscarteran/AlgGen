!--------------------------------------------------------------------------------------------------!
! Función conv1d.f08 [Nava-Flores, 2020]
!
! Útil para realizar la convolución discreta en 1D.
!
! Parámetros de entrada:
!  - Arreglo por convolucionar (entrada): x <- Arreglo 1D
!  - Kernel convolutivo: h <------------------ Arreglo 1D
!  - Marcador para solicitar la salida con el mismo numero de elementos que el arreglo de entrada 
!    (opcional): ind <------------------------ Cadena de caracteres
!                ind = 'same' -> La salida tiene la misma longitud que la entrada
!
! Últimas modificaciones: Nava-Flores, 23/Ene/2021
!--------------------------------------------------------------------------------------------------!
function conv1d(x, h, ind)

   use NumKind

   implicit none
   character(len=4), optional, intent(in):: ind
   integer(il):: i, j, N1, N2, N3
   real(dp), dimension(:), intent(in):: x, h
   real(dp), allocatable, dimension(:):: conv1ds, conv1d

   ! Dimensiones del arreglo de entrada:
   N1=size(x)

   ! Dimensiones del kernel:
   N2=size(h)

   ! Dimensiones del arreglo de salida:
   N3=N1+N2-1

   ! Arreglo de salida:
   allocate(conv1ds(N3))

   conv1ds=0._dp
   
   do i=1,N1
      do j=0,N2-1
         conv1ds(i+j)=conv1ds(i+j) + x(i)*h(j+1)
      end do
   end do

   if (N2 > N1) then
      i=N1
      N1=N2
      N2=i
   end if

   if (present(ind)) then
      if (ind == 'same' .or. ind == 'SAME' .or. ind == 'Same') then
         allocate(conv1d(N1))
         conv1d = conv1ds(floor(real(N2, dp)/2)+1:floor(real(N2, dp)/2)+N1)
      else
         print '(a)', 'Parámetro "ind" erróneo'
         stop
      end if
   else
      allocate(conv1d(N3))
      conv1d = conv1ds
   end if

end function conv1d
!--------------------------------------------------------------------------------------------------!