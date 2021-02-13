!--------------------------------------------------------------------------------------------------!
! Módulo NumKind [Nava-Flores, 2017]
! - Basado en el módulo NumericKinds [Chirila, D. B., & Lohmann, G. (2015)].
!
! Calcula y establece los números "KIND" y sus "alias", para los siguientes tipos de variables:
!  - Integer:
!     > Simple (IS)
!     > Long (IL)
!     > Long Long (ILL)
!
!  - Real:
!     > Single Precision (SP)
!     > Double Precision (DP)
!     > Quadruple Precision (QP)
!
! El objetivo de este módulo es definir la precisión de este tipo de variables, independientemente 
! del compilador utilizado.
!
! Una vez calculados los números kind, se pueden utilizar a través de sus "alias" 
! (IS, IL, ILL, SP, DP y QP)
!
! Última modificación: Nava-Flores, 26/Oct/2020
!--------------------------------------------------------------------------------------------------!
module NumKind

   implicit none

   !Números KIND para diferentes tipos de enteros:
   integer, parameter:: IS = selected_int_kind(4), &
                        IL = selected_int_kind(9), &
                        ILL= selected_int_kind(18)

   !Números KIND para diferentes tipos de reales:
   integer, parameter:: SP = selected_real_kind(6,37), &
                        DP = selected_real_kind(15,307), &
                        QP = selected_real_kind(33,4931)

end module NumKind
!--------------------------------------------------------------------------------------------------!