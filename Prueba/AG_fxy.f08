!--------------------------------------------------------------------------------------------------!
! UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
! FACULTAD DE INGENIERÍA
! DIVISIÓN DE INGENIERÍA EN CIENCIAS DE LA TIERRA
! DEPARTAMENTO DE GEOFÍSICA
!
! Asignatura: Inversión de Datos Geofísicos
!   Profesor: Dr. Mauricio Nava Flores
!--------------------------------------------------------------------------------------------------!
! Programa en Fortran 2008 que implementa un algoritmo genético para determinar el mínimo global de
! una función 2D: f(x,y).
!
! La estructura del programa se basa en el algoritmo propuesto por:
! 
!     Monte Carlo Methods for strongly non‐linear geophysical optimization problems. Geophysical 
!     Gallagher, K., Sambridge, M., & Drijkoningen, G. (1991). Genetic algorithms: An evolution from
!     Research Letters, 18(12), 2177-2180.
! 
!------------------------------------------------
! Árbol de dependencias del programa AG_fxy.f08:
! 
!  AG_fxy -> NumKind
!         -> dec2bin -> NumKind
!         -> bin2dec -> NumKind
!         -> mean ----> NumKind
!         -> init_random_seed
!------------------------------------------------
!
! Últimas modificaciones: Nava-Flores, 20/Ene/2021
!--------------------------------------------------------------------------------------------------!
program AG_fxy
   
   use NumKind
   use modRW
   
   implicit none
   
   integer(il):: i, j, Q, M, N, NGen, NR, NBits, iGen, NApca, unit1, Nava, unit2, unit3, it
   
   real(dp):: t1, xmin, xmax, ymin, ymax, Pc, Pm, Tol, phim, phii, b, a, &
              rand1, NAPc, NAPm, tempi, tempj, t2
   
   character(len=100):: DirRes, filename
   
   character(len=:), allocatable:: tempB1, tempB2, tempB1H, tempB2H
   
   character(len=200), allocatable:: chBits(:)
   
   integer(il), allocatable, dimension(:):: mMax, Bits, Ind, tempInd1, tempInd2, tempiGeni,       &
                                            tempiGenj, cont1
   
   integer(il), allocatable, dimension(:,:):: IndGen, IndGenP, IndGenRep
   
   real(dp), allocatable, dimension(:):: pMin, pMax, dm, phi, randM, Pr, randQ, NAPr, GenOpt,     &
                                         tempGeni, tempGenj, DatX, DatY, tconv, f, Zobs, Zest,  &
                                         Zopti, Zopt
   
   real(dp), allocatable, dimension(:,:):: EspM, Gen, GenP, GenRep, Dat
   
   real(dp), parameter:: pi=4._dp*atan(1._dp)
   

   !----- Interfaz explícita para incluir los subprogramas requeridos:
   interface
      function dec2bin(d,n)
         use NumKind
         integer(il), intent(in):: d
         integer(il), intent(in), optional:: n
         character(len=:), allocatable:: dec2bin
      end function dec2bin

      function bin2dec(cadenabin)
         use NumKind
         character(len=*), intent(in):: cadenabin
         integer(il):: bin2dec
      end function bin2dec

      function mean(x)
         use NumKind
         real(dp),allocatable,intent(in):: x(:)
         real(dp):: mean
      end function mean

      subroutine init_random_seed()
         use iso_fortran_env, only: int64
      end subroutine init_random_seed
   end interface
   !-----
   
   !Lectura de datos:
   Dat = readMat('AnomG_1.dat')
   DatX = Dat(:,1)
   Zobs = Dat(:,2)
   N = size(Dat,1) 
   
   !Lectura de Datos 2D
   !DatY = Dat(:,2)
   !Zobs = Dat(:,3)

   !----- Función por optimizar:
   !f(x,y) = x**2 + y**2
   !----- Anomalía Gravimetrica: 12 Parametros
   !Zobs = CilindroG(x,a,xc,zc,rho)+CilindroG(x,a,xc,zc,Rho) + EsferaG(x,a,xc,zc,Rho)
   !----- Anomalia Magnetica: 
   !Zobs = EsferaM(x,y,x,y,z,a,mi,md,mg)
   !----- Sismica
   !Zobs = conv1d(tconv,Ricker,'same')
   !----- Eléctrica
   !Zobs = SEV(ab2,m,3)
   
   !-----

   ! Limpieza de terminal:
   call system('clear')

   ! Inicia cronómetro:
   call cpu_time(t1)

   ! Directorio de resultados:
   DirRes='/home/oscar/Escritorio/AlgGen/Prueba/DirResl/'

   !-----------------------------------------------------------------------------------------------!
   !                                    FUNCIÓN POR OPTIMIZAR
   !-----------------------------------------------------------------------------------------------!
   !xmin=-5.12_dp
   !xmax= 5.12_dp
   !ymin=-5.12_dp
   !ymax= 5.12_dp
   
   !Definir espacio de busqueda con vectores xmin y x max
   
   !xmin = minval(x)
   !xmax = maxval(x)
   !ymin = minval(y)
   !ymax = maxval(y)

   ! Mínimo absoluto de la función f(x,y):
   !Zmin=0.0_dp
   !-----------------------------------------------------------------------------------------------!

   !-----------------------------------------------------------------------------------------------!
   !                      PARÁMETROS CONTROLADORES DEL ALGORITMO GENÉTICO
   !-----------------------------------------------------------------------------------------------!
   ! Parámetros Generales:
   Q = 100_il                              ! Tamaño de población
   M = 12_il                               ! Número de parámetros, dependiendo de la inversión
   Pc = 1._dp                            ! Probabilidad de cruza (Cte: 90%)
   Pm = 0.1_dp                           ! Probabilidad de mutación (Cte: 5%)
   NGen = 10000_il                        ! Número máximo de generaciones
   Tol = 0.01_dp                         ! Tolerancia (c/R al error)
   !-----------------------------------------------------------------------------------------------!
   
   
   ! Inicialización del generador de números aleatorios:
   call init_random_seed()
    
   allocate(pMin(M), pMax(M), dm(M), mMax(M), chBits(M), Bits(M))
   ! Espacio discreto de modelos:
      !Gravimetria
   !pMin = [0.01_dp, 0.170_dp, 0.020_dp, -2000._dp, &
   !	   0.02_dp, 0.360_dp, 0.020_dp, 600._dp, &
   !	   0.05_dp, 0.710_dp, 0.020_dp, 500._dp]!posiciones de los minimos
   
   !pMax = [0.06_dp, 0.200_dp, 0.080_dp, -1500._dp, &
   !	 0.05_dp, 0.400_dp, 0.080_dp, 3500._dp, &
   !	  0.07_dp, 0.850_dp, 0.080_dp, 3500._dp]!posiciones de los maximos
   !pMin = [0.007_dp, 0.170_dp, 0.060_dp, -15000._dp, &
   !      0.005_dp, 0.340_dp, 0.060_dp, 10000._dp, &
   !      0.010_dp, 0.750_dp, 0.040_dp, 25000._dp]!posiciones de los minimos
   
   !pMax = [0.01_dp, 0.200_dp, 0.080_dp, 10000._dp, &
   !    0.010_dp, 0.400_dp, 0.080_dp, 15000._dp, &
   !     0.020_dp, 0.775_dp, 0.080_dp, 35000._dp]!posiciones de los maximos
   pMin = [0.01_dp, 0.180_dp, 0.012_dp, -3500._dp, &
         0.02_dp, 0.370_dp, 0.020_dp, 600._dp, &
         0.05_dp, 0.750_dp, 0.050_dp, 500._dp]!posiciones de los minimos
   
   pMax = [0.06_dp, 0.200_dp, 0.080_dp, -1500._dp, &
       0.05_dp, 0.400_dp, 0.080_dp, 3500._dp, &
        0.07_dp, 0.780_dp, 0.080_dp, 2000._dp]!posiciones de los maximos
   
   dm = 30.0_dp
   mMax = nint((pMax-pMin)/dm+1)
   
   !suma total de los posibles parametros
   allocate(EspM(maxval(mMax),M))	  
   EspM=0._dp
   do i=1,M
      do j=1,mMax(i)
         EspM(j,i)=pmin(i)+(j-1)*dm(i)
      end do
   end do
   
   ! Número de bits para la codificación binaria:
   write(chBits, '(b0)') mMax
   Bits = len_trim(chBits) !Porque marcas error?
   NBits = sum(Bits) 
   !-----------------------------------------------------------------------------------------------!
   !  INICIO DEL ALGORITMO:
   !-----------------------------------------------------------------------------------------------!
   allocate(Zest(N))
   allocate(Gen(Q,M), GenOpt(M), IndGen(Q,M), phi(Q)) !Porque marcas error
   allocate(randM(M), Ind(M))
   allocate(Pr(Q), GenP(Q,M), IndGenP(Q,M), randQ(Q))
   allocate(GenRep(Q,M), IndGenRep(Q,M), tempInd1(M), tempInd2(M))
   allocate(character(NBits) :: tempB1, tempB2, tempB1H, tempB2H)
   allocate(tempGeni(M), tempGenj(M), tempiGeni(M), tempiGenj(M))
   allocate(cont1(Q))
   
   !open(newunit=unit2, file=trim(DirRes)//'Resultados.dat', status='replace', action='write')
   
   Gen=0._dp                ! Arreglo 2D para almacenar los "Q" modelos
   GenOpt=0._dp             ! Arreglo 1D para almacenar el modelo óptimo
   IndGen=0_il              ! Arreglo 2D para almacenar los índices de los modelos
   phi=0._dp                ! Arreglo 1D para almacenar el "error" de los modelos
   phim=1000._dp            ! Valor inicial de error mínimo global (se actualizará solo)
   phii=1000._dp            ! Valor inicial de error mínimo de c/generación (se actualizará solo)
   !-----------------------------------------------------------------------------------------------!
   
   


   

  

   ! Almacenamiento de la curva de convergencia en archivo ASCII:
   !WRITE (filename, fmt='(a,I0,a)') &
   !'/home/alejandro/Desktop/IDGB/P_Final/DirResl/CurvaConv',Nava,'.dat'
   !OPEN (UNIT = unit1, FILE=filename, STATUS = 'REPLACE', action='write')
   !-----------------------------------------------------------------------------------------------!

   !-----------------------------------------------------------------------------------------------!
   !  CICLO SOBRE EL NÚMERO MÁXIMO DE GENERACIONES (O HASTA ALCANZAR LA TOLERANCIA ESTABLECIDA)
   !-----------------------------------------------------------------------------------------------!
   iGen=1
   
   do while (iGen <= NGen)
      !Respuesta de la población al medio:
      do i=1,Q
         if (iGen == 1) then
            ! Índices aleatorios para cada parámetro:
            call random_number(randM)
            Ind = int(1 + (mMax-1)*randM)
            IndGen(i,:) = Ind
         else
            ! Índices heredados de la generación anterior:
            Ind = IndGen(i,:)
         end if

         ! Genes de los individuos:
         do j=1,M
            Gen(i,j) = EspM(Ind(j),j)
         end do

         ! Respuesta de los individuos al medio:
         !Medio Gravimetria
         do it =1,N
         Zest(it) = CilindroG(DatX(it), Gen(i,1), Gen(i,2), Gen(i,3), Gen(1,4)) + &
         	CilindroG(DatX(it), Gen(i,5), Gen(i,6), Gen(i,7), Gen(i,8)) + &
         	EsferaG(DatX(it), Gen(i,9), Gen(i,10), Gen(i,11), Gen(i,12))
         end do
         !Zest = f(Gen(i,1), Gen(i,2))
	
	!= cilindro()+cilindro()+esfera()
	!gen->(i,12)
	
	!vector zestimada
         ! Cálculo del desajuste:
         phi(i) = norm2(Zobs - Zest)

         ! Búsqueda del individuo más apto de la generación actual:
         if (iGen == 1) then
            if (phi(i) < phii) then
               phii = phi(i)
               phim = phi(i)
               Zopt = Zest
               Zopti = Zest
               GenOpt = Gen(i,:)
            end if
         elseif (iGen > 1) then
            if (phi(i) < phii) then
               phii = phi(i)
               Zopti = Zest
               ! *** Si el individuo más apto de la generación actual está más adaptado al medio que
               !     el  óptimo hasta el momento, se tendrá un nuevo óptimo global:
               if (phi(i) < phim) then
                  phim = phi(i)
                  Zopt = Zest
                  GenOpt = Gen(i,:)
               end if
            end if
         end if
      end do

      ! Convergencia (iGen vs phi):
      !write(unit1,*) iGen, phim
      !--------------------------------------------------------------------------------------------!

      !Criterio de detención del algoritmo:
      if (phim <= Tol) exit

      !--------------------------------------------------------------------------------------------!
      !  REPRODUCCIÓN:        SELECCIÓN DE INDIVIDUOS POR RULETA
      !--------------------------------------------------------------------------------------------!
      ! Probabilidad de reproducción en función del error (phi):
      Pr=0._dp
      GenP=0._dp
      IndGenP=0_il
      b = -1._dp/(Q*(maxval(phi)-mean(phi)))
      a = -b*maxval(phi)
      do i=1,Q
         Pr(i) = a + b*phi(i)
      end do

      ! Probabilidades acumuladas:
      do i=1,Q
         if (i == 1) then
            Pr(i) = Pr(i)
         else
            Pr(i) = Pr(i) + Pr(i-1)
         end if
      end do

      ! Selección de individuos para la reproducción por ruleta:
      call random_number(randQ)
      NAPr = randQ
      do i=1,Q
         do j=1,Q
            if (j == 1) then
               if (NAPr(i) < Pr(j)) then
                  GenP(i,:) = Gen(j,:)
                  IndGenP(i,:) = IndGen(j,:)
               end if
            elseif (j > 1 .and. j < Q) then
               if (NAPr(i) >= Pr(j-1) .and. NAPr(i) < Pr(j)) then
                  GenP(i,:) = Gen(j,:)
                  IndGenP(i,:) = IndGen(j,:)
               end if
            elseif (j == Q) then
               if (NAPr(i) >= Pr(j-1)) then
                  GenP(i,:) = Gen(j,:)
                  IndGenP(i,:) = IndGen(j,:)
               end if
            end if
         end do
      end do
      !--------------------------------------------------------------------------------------------!

      !--------------------------------------------------------------------------------------------!
      !  REPRODUCCIÓN:        CRUZA CONTROLADA POR Pc Y MUTACIÓN CONTROLADA POR Pm
      !--------------------------------------------------------------------------------------------!
      GenRep=0._dp
      IndGenRep=0_il
      tempInd1=0_il
      tempInd2=0_il

      do i=1,Q,2
         call random_number(rand1)
         NAPc = rand1

         ! Si NAPc <= Pc, los individuos se cruzan:
         if (NAPc <= Pc) then
            do j=1,M
               tempInd1(j) = IndGenP(i,j)
               tempInd2(j) = IndGenP(i+1,j)
            end do

            ! Codificación binaria:
            tempB1(1:Bits(1)) = dec2bin(tempInd1(1),Bits(1))
            tempB1(Bits(1)+1:NBits) = dec2bin(tempInd1(2),Bits(2))
            tempB2(1:Bits(1)) = dec2bin(tempInd2(1),Bits(1))
            tempB2(Bits(1)+1:NBits) = dec2bin(tempInd2(2),Bits(2))

            !----- Cruza en punto aleatorio:
            call random_number(rand1)
            NApca = nint(1 + (NBits-1)*rand1)
            ! Descendiente 1:
            tempB1H(1:NApca) = tempB1(1:NApca)
            tempB1H(NApca+1:NBits) = tempB2(NApca+1:NBits)
            ! Descendiente 2:
            tempB2H(1:NApca) = tempB2(1:NApca)
            tempB2H(NApca+1:NBits) = tempB1(NApca+1:NBits)
            !-----

            !----- Mutación controlada por Pm (Cte):
            call random_number(rand1)
            NAPm = rand1
            if (NAPm <= Pm) then
               ! Mutación en bit aleatorio - Descendiente 1:
               call random_number(rand1)
               NApca = nint(1 + (NBits-1)*rand1)
               if (tempB1H(NApca:NApca) == '1') then
                  tempB1H(NApca:NApca) = '0'
               else
                  tempB1H(NApca:NApca) = '1'
               end if
            end if

            call random_number(rand1)
            NAPm = rand1
            if (NAPm < Pm) then
               ! Mutación en bit aleatorio - Descendiente 2:
               call random_number(rand1)
               NApca = nint(1 + (NBits-1)*rand1)
               if (tempB2H(NApca:NApca) == '1') then
                  tempB2H(NApca:NApca) = '0'
               else
                  tempB2H(NApca:NApca) = '1'
               end if
            end if
            !-----

            !----- Decodificación y chequeo de límites:
            do j=1,M
               tempInd1(j) = bin2dec(tempB1H((j-1)*NBits/M+1:j*NBits/M))
               tempInd2(j) = bin2dec(tempB2H((j-1)*NBits/M+1:j*NBits/M))
               if (tempInd1(j) < 1 .or. tempInd1(j) > mMax(j)) then
                  call random_number(rand1)
                  tempInd1(j) = nint(1 + (mMax(j)-1)*rand1)
               end if
               if (tempInd2(j) < 1 .or. tempInd2(j) > mMax(j)) then
                  call random_number(rand1)
                  tempInd2(j) = nint(1 + (mMax(j)-1)*rand1)
               end if
            end do
            IndGenRep(i,:) = tempInd1
            IndGenRep(i+1,:) = tempInd2
            !-----

         else   
            ! Si NAPc > Pc, los individuos se clonan (no hay cruza):
            do j=1,M
               IndGenRep(i,j) = IndGenP(i,j)
               IndGenRep(i+1,j) = IndGenP(i+1,j)
            end do
         end if
      end do
      !--------------------------------------------------------------------------------------------!

      !--------------------------------------------------------------------------------------------!
      !  SELECCIÓN CON ELITISMO Y REEMPLAZO DE INDIVIDUOS IDÉNTICOS
      !--------------------------------------------------------------------------------------------!
      ! Respuesta de la nueva generación al medio:
      do i=1,Q
         ! Índices de la nueva generación:
         Ind = IndGenRep(i,:)

         ! Genes de los individuos:
         do j=1,M
            Gen(i,j) = EspM(Ind(j),j)
         end do

         ! Respuesta de los individuos al medio:
         !Zest = f(Gen(i,1),Gen(i,2))
         !Medio Gravimétrico:
         do it =1,N
         Zest(it) = CilindroG(DatX(it), Gen(i,1), Gen(i,2), Gen(i,3), Gen(1,4)) + &
         	CilindroG(DatX(it), Gen(i,5), Gen(i,6), Gen(i,7), Gen(i,8)) + &
         	EsferaG(DatX(it), Gen(i,9), Gen(i,10), Gen(i,11), Gen(i,12))
         end do

         ! Cálculo del desajuste:
         phi(i) = norm2(Zobs - Zest)
      end do

      ! Orden de la nueva población c/r al desajuste (menor a mayor):
      tempGeni=0._dp
      tempGenj=0._dp
      tempiGeni=0_il
      tempiGenj=0_il
      do i=1,Q
         do j=1,Q
            tempi = phi(i);   tempGeni = Gen(i,:);   tempiGeni = IndGenRep(i,:)
            tempj = phi(j);   tempGenj = Gen(j,:);   tempiGenj = IndGenRep(j,:)
            if (tempi < tempj) then
               phi(i) = tempj;  Gen(i,:) = tempGenj;   IndGenRep(i,:) = tempiGenj
               phi(j) = tempi;  Gen(j,:) = tempGeni;   IndGenRep(j,:) = tempiGeni
            end if
         end do
      end do

      ! El individuo más adaptado pasa sin alterarse a la nueva generación (elitismo):
      IndGen=0
      IndGen(1,:) = IndGenRep(1,:)

      ! Conteo del número de individuos idénticos:
      cont1=0
      do j=1,Q
         do i=j+1,Q
            if (sum(abs(IndGenRep(j,:)-IndGenRep(i,:))) == 0) then
               cont1(j) = cont1(j) + 1
            end if
         end do
      end do

      ! Nueva generación permitiendo hasta NR modelos idénticos repetidos:
      NR=Q/4
      do i=1,Q
         if (cont1(i) <= NR) then
            IndGen(i,:) = IndGenRep(i,:)
         elseif (cont1(i) > NR) then
            call random_number(randM)
            IndGen(i,:)=int(1+(mMax-1)*randM)
         end if
      end do

      ! Despliegue de resultados parciales en terminal:
      write(*,*)
      write(*,100) iGen, phim

      ! Fin de iteración actual e inicio de la siguiente:
      phii=1000._dp
      iGen=iGen+1
   end do

   !close(unit1)

   ! Se detiene el cronómetro:
   call cpu_time(t2)
   
   !Se almacenan los vectores de parámetros y anomalía invertidos:
   open(NewUnit=unit2,file=trim(DirRes)//'ParamOpts.dat',status='replace',action='write')
   do i=1,M
      write(unit2,*) GenOpt(i)
   end do
   close(unit2)

   open(NewUnit=unit3,file=trim(DirRes)//'AnomOpt.dat',status='replace',action='write')
   do i=1,N
      write(unit3,*) DatX(i), Zopt(i)
   end do
   close(unit3)

   ! Despliegue de resultados generales de la inversión:
   write(*,110)
   write(*,120) M
   write(*,130) iGen-1
   write(*,140) Tol
   write(*,150) phim
   !write(*,160) Zmin
   !write(*,170) Zopt
   !write(*,180) GenOpt(1),GenOpt(2)
   write(*,190) Q
   write(*,200) Pc
   write(*,210) Pm
   write(*,220) NBits
   write(*,230) t2-t1

   100 format('----------------------------------------------------------------------',/,&
              '            RESULTADOS PARCIALES EN LA GENERACIÓN ',i0,'              ',/,&
              ' Desajuste mínimo de la generación actual: ',E15.8,'                  ',/,&
              '----------------------------------------------------------------------')

   110 format(/,'----------------------------------------------------------------------'/&
                '         RESULTADOS DE LA INVERSIÓN POR ALGORITMOS GENÉTICOS          '/&
                '----------------------------------------------------------------------')
   120 format('       Número de parámetros invertidos: 'i0)
   130 format('                Número de generaciones: 'i0)
   140 format('                            Tolerancia: 'f7.4)
   150 format('            Desajuste mínimo alcanzado: 'f9.6)
   160 format('                Mínimo global conocido: 'f9.6)
   170 format('               Mínimo global calculado: 'f9.6)
   180 format('      Coordenadas del óptimo calculado: ('f7.4','f7.4')'/)
   190 format(' Número de genomas en la población (Q): 'i0)
   200 format('            Probabilidad de cruza (Pc): 'f4.2)
   210 format('         Probabilidad de mutación (Pm): 'f4.2)
   220 format('         Número de Bits para codificar: 'i0/)
   230 format('                     Tiempo de cómputo: 'f0.4' [seg]'/,&
              '------------------------------------------------------------------------'/)
   !write(unit2,*) Nava, iGen-1, Zopt, GenOpt(1), GenOpt(2), t2-t1
!end do
!close(unit2)

end program AG_fxy
!--------------------------------------------------------------------------------------------------!
