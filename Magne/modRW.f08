!--------------------------------------------------------------------------------------------------!
! UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
! FACULTAD DE INGENIERÍA
! DIVISIÓN DE INGENIERÍA EN CIENCIAS DE LA TIERRA
! DEPARTAMENTO DE GEOFÍSICA
!
! Asignatura: Inversión de Datos Geofísicos
!   Profesor: Dr. Mauricio Nava Flores
!
! Módulo modRW empleado para englobar funciones de lectura y escritura de archivos ASCII (1D y 2D).
!
! El módulo ModRW contiene las siguientes funciones:
!  1. readVect. Sirve para leer datos de archivos ASCII de una columna sin encabezados.
!  2. readMat. Sirve para leer datos de archivos ASCII de más de una columna sin encabezados.
! 
! El módulo ModRW contiene las siguientes subrutinas:
!  1. writeVect. Sirve para escribir datos de un vector (arreglo 1D) en archivo ASCII sin 
!                encabezados.
!  2. writeMat. Sirve para escribir datos de una matriz (arreglo 2D) en archivo ASCII sin 
!               encabezados.
! 
! Fecha de entrega: 11 de diciembre de 2020
! 
! Autores:
!           Hernández Terán Oscar 
!           Romero Reyes Alejandro
!           Martinez Reyes Javier
!--------------------------------------------------------------------------------------------------!
module modRW

   use NumKind
   contains

! --------------------------------------------------------------------------------------------------
! Función para la lectura de un vector en un archivo tipo ASCII
! --------------------------------------------------------------------------------------------------
   function readVect(Arch) result(x)
      
      implicit none
      character(len=*):: Arch
      integer(il):: i, N, UnLec, stat
      real(dp), allocatable:: x(:)

      ! Apertura de unidad de lectura:
      open(NewUnit=UnLec, file=Arch, status='old', action='read')

      ! Determinación del número de datos (renglones o filas) del archivo:
      N=0
      do 
         read(UnLec,*,iostat=stat)

         ! El ciclo se detiene al momento en que sea imposible seguir leyendo:
         if (stat /= 0) exit
         N=1+N
      end do

      ! El puntero de lectura se regresa al inicio del archivo:
      rewind(UnLec)

      ! Asignación de memoria al arreglo x:
      allocate(x(N))

      ! Lectura de datos en el arreglo x:
      do i=1,N
         read(UnLec,*) x(i)
      end do

      ! Se cierra la unidad de lectura:
      close(UnLec)

      ! Se indica el número de datos leídos:

      print '(a)','---------------------------------------'
      print '(a,i0,/)', 'Número de datos leídos: ', N
      print '(a)','---------------------------------------'

      

   end function readVect

! --------------------------------------------------------------------------------------------------
! Función para la lectura de una matriz en un archivo tipo ASCII
! --------------------------------------------------------------------------------------------------
  function readMat(Arch) result(x)
      
      implicit none
      logical:: delimitador=.false.
      character(len=*):: Arch
      character(len=1000000):: buffer
      integer:: i, j, M, N, UnLec, stat
      double precision, allocatable:: x(:,:)

      ! Apertura de unidad de lectura:
      open(NewUnit=UnLec, file=Arch, status='old', action='read')

      ! Determinación del número de renglones o filas del archivo:
      M=0
      do 
         read(UnLec,*,iostat=stat)

         ! El ciclo se detiene al momento en que sea imposible seguir leyendo:
         if (stat /= 0) exit
         M=1+M
      end do
      print '(a)','---------------------------------------'
      print '(a,i0)', 'Número de renglones: ', M

      ! El puntero de lectura se regresa al inicio del archivo:
      rewind(UnLec)

      !----- Determinación del número de columnas del archivo (con una fila de datos):
      read(UnLec,'(a)',iostat=stat) buffer
      if (stat /= 0) then
         print '(a)', 'Error de lectrura!!!'
      end if

      buffer=adjustl(buffer)

      N=1
      do i=1,len_trim(buffer)
         if (ichar(buffer(i:i))==09 .or. ichar(buffer(i:i))==44 .or. ichar(buffer(i:i))==32) then
            if (delimitador .eqv. .false.) then
               N=N+1
               delimitador=.true.
            end if
         else
            delimitador=.false.
         end if
      end do

      i=len_trim(buffer)
      if (ichar(buffer(i:i))==44 .or. ichar(buffer(i:i))==09) N=N-1
      print '(a,i0,/)', 'Número de columnas: ', N
      print '(a)','---------------------------------------'
      !-----

      ! El puntero de lectura se regresa al inicio del archivo:
      rewind(UnLec)

      ! Asignación de memoria al arreglo x:
      allocate(x(M,N))

      ! Lectura de datos en el arreglo x:
      do i=1,M
         read(UnLec,*) (x(i,j), j=1,N)
      end do

      ! Se cierra la unidad de lectura:
      close(UnLec)

   end function readMat

! --------------------------------------------------------------------------------------------------
! Función para la escritura de un vector en un archivo tipo ASCII
! --------------------------------------------------------------------------------------------------
   subroutine writeVect(Arch, Vector)

      implicit none
      character(len=*):: Arch
      integer(il):: UnEsc, i, N
      real(dp):: Vector(:)

      ! Definición del número de elementos del vector "Vector":
      N=size(Vector,1)

      ! Apertura de la unidad de escritura de datos en formato ASCII:
      open(NewUnit=UnEsc, file=Arch, status='replace', action='write')
      do i=1,N
         write(UnEsc,'(f0.12)') Vector(i)
      end do
      close(UnEsc)

      print '(a,i0,/)', 'Número de datos escritos: ', N

   end subroutine writeVect

! --------------------------------------------------------------------------------------------------
! Función para la escritura de una matriz en un archivo tipo ASCII
! --------------------------------------------------------------------------------------------------

   subroutine writeMat(Arch, Matriz)

      implicit none
      character(len=*):: Arch
      integer(il):: UnEsc, i, j, N, M
      real(dp):: Matriz(:,:)

      ! Definición del número de elementos de la matriz "Matriz":
      M=size(Matriz,1)
      N=size(Matriz,2)

      ! Apertura de la unidad de escritura de datos en formato ASCII:
      open(NewUnit=UnEsc, file=Arch, status='replace', action='write')
      do i=1,M
         write(UnEsc,*) (Matriz(i,j), j=1,N)
      end do
      close(UnEsc)

   end subroutine writeMat

! --------------------------------------------------------------------------------------------------

! --------------------------------------------------------------------------------------------------
! FUNCIONES PARA INVERTIR DATOS GEOFISICOS
! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------
! CILINDRO GRAVIMETRÍA
! --------------------------------------------------------------------------------------------------
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

   EsferaG=(-k*(0-zc)/(R**3))*si2mg*km2m
   
end function EsferaG
!--------------------------------------------------------------------------------------------------!

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
function EsferaM(x0, y0, xc, yc, zc, a, mi, md, mg)

   use NumKind

   implicit none
   real(dp), intent(in):: x0, y0, xc, yc, zc, a, mi, md, mg
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
   rz=-zc
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

!--------------------------------------------------------------------------------------------------
!Se genera el vector de reflectividades
!--------------------------------------------------------------------------------------------------
function ModR(Amp) result(y)
   use NumKind

   implicit none
   integer(il):: i, j, Ind(7)
   real(dp):: Amp(7), y(1001)

   ! Ejemplo: Amp=[-3._dp, -2._dp, ..., 1.75_dp]

   Ind=[101, 201, 301, 501, 601, 701, 901]
   
   y=0._dp
   do i=1,1001
      do j=1,7
         if (i == Ind(j)) y(i) = Amp(j)
      end do
   end do
end function ModR

!-----------------------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------------------!
! Función que calcula el problema directo para un SEV adquirido en arreglo electródico Schlumberger
! en un punto de apertura (ab/2), considerando un semiespacio estratificado en capas infinitamente
! extendidas.
! 
!  - El cálculo del problema directo se hace a través de una convolución, con base en la teoría de 
!    filtrado lineal (Ghosh, 1971), empleando los coeficientes y estructura reportados por Ekinci y 
!    y Demirci (2008).
! 
! Parámetros de entrada:
!  - ab2: Punto de abertura ab/2 del arreglo electródico Schlumberger
!  -   m: Vector de resistividades (rho) y espesores (h) del medio
!           m = rho_1, rho_2, ..., rho_nc, h_1, h_2, ..., h_nc-1
!  -  nc: Número de capas del medio
!
! Referencias:
! 
!     Ghosh, D. P. (1971). The application of linear filter theory to the direct interpretation of 
!     geoelectrical resistivity sounding measurements. Geophysical prospecting, 19(2), 192-217.
!
!     Ekinci, Y. L., & Demirci, A. (2008). A damped least-squares inversion program for the 
!     interpretation of Schlumberger sounding curves. Journal of Applied Sciences, 8(22), 4070-4078.
!
! Últimas modificaciones: Nava-Flores, 23/Ene/2021
!--------------------------------------------------------------------------------------------------!
function SEV(ab2, m, nc)

   use NumKind

   implicit none
   integer(il), intent(in):: nc
   integer(il):: i, j, q
   real(dp), intent(in):: ab2
   real(dp), intent(in):: m(:)
   real(dp):: rho(nc), h(nc-1), c(25), R, mm
   real(dp):: SEV

   rho=m(1:nc)
   h=m(nc+1:2*nc-1)

   ! Orden del filtro:
   q=13_il

   ! Parámetro del filtro relacionado con el muestreo:
   mm=4.438_dp

   ! Coeficientes del filtro:
   c=[105._dp, 0._dp, -262._dp, 0._dp, 416._dp, 0._dp, -746._dp, 0._dp, 1605._dp, 0._dp, -4390._dp,&
      0._dp, 13396._dp, 0._dp, -27841._dp, 0._dp, 16448._dp, 0._dp, 8183._dp, 0._dp, 2525._dp,     &
      0._dp, 336._dp, 0._dp, 225._dp]
   c=c/10000._dp
   
   SEV=0._dp
   do i=1,2*q-1
      R=rho(nc)
      do j=nc,2,-1
         R=(R+rho(j-1)*tanh(h(j-1)/(ab2*exp(((i-1)/2._dp-10._dp)*log(10._dp)/mm))))/ &
           (1._dp+R/rho(j-1)*tanh(h(j-1)/(ab2*exp(((i-1)/2._dp-10._dp)*log(10._dp)/mm))))
      end do
      SEV=SEV+c(i)*R
   end do

end function SEV
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
! Subrutina que genera una señal con la forma del pulso Ricker.
! 
! Parámetro de entrada:
!  - fm: Frecuencia central
! 
! Parámetro de salida:
!  - t: Vector de tiempo
!  - f: Pulso Ricker 
!
! Últimas modificaciones: Nava-Flores, 23/Ene/2021
!--------------------------------------------------------------------------------------------------!
subroutine Ricker(fm, t, f)

   use NumKind

   implicit none
   integer(il):: i, N
   real(dp), intent(in):: fm
   real(dp):: tini, tfin, dt, t0, a
   real(dp), allocatable:: t(:), f(:)
   real(dp), parameter:: pi=4._dp*atan(1._dp)

   dt=1._dp/(20._dp*fm)
   a=(pi*fm)**2

   tini=0._dp
   tfin=2.5_dp/fm
   t0=1.25_dp/fm

   N=nint((tfin-tini)/dt + 1, il)

   allocate(t(N), f(N))	
   t=[(tini + (i-1)*dt, i=1, N)]

   f=-(2._dp*a*((t-t0)**2)-1._dp)*(exp(-a*((t-t0)**2)))

end subroutine Ricker
!--------------------------------------------------------------------------------------------------!

end module modRW

!--------------------------------------------------------------------------------------------------!
