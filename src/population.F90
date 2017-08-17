#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: deb_population
!
! !INTERFACE:
module deb_population
!
! !DESCRIPTION:
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_population
      ! Variable identifiers
      type (type_bottom_state_variable_id), allocatable :: id_NV(:)
      type (type_bottom_state_variable_id), allocatable :: id_NE(:)
      type (type_bottom_state_variable_id), allocatable :: id_NE_H(:)
      type (type_bottom_state_variable_id), allocatable :: id_NQ(:)
      type (type_bottom_state_variable_id), allocatable :: id_NH(:)

      type (type_bottom_state_variable_id) :: id_waste

      type (type_horizontal_diagnostic_variable_id) :: id_R, id_N_e, id_N_j, id_N_a

      type (type_dependency_id) :: id_temp
      type (type_dependency_id) :: id_salt

      ! Number of size classes and prey
      integer :: nclass

      ! Size class characteristics
      real(rk),allocatable :: logV_center(:)     ! log structural volume
      real(rk),allocatable :: V_center(:)        ! structural volume
      real(rk),allocatable :: delta_V(:)         ! structural length difference between consecutive size classes
      real(rk),allocatable :: L_center(:)        ! structural length

      ! Size-class-independent parameters
      real(rk) :: p_Am
      real(rk) :: p_T
      real(rk) :: p_M
      real(rk) :: E_G
      real(rk) :: E_Hb
      real(rk) :: E_Hj
      real(rk) :: E_Hp
      real(rk) :: v
      real(rk) :: k_J
      real(rk) :: kap
      real(rk) :: kap_R
      real(rk) :: h_a
      real(rk) :: s_G

      real(rk) :: E_0
      real(rk) :: L_T
      real(rk) :: E_m
      real(rk) :: s_M
      real(rk) :: L_b, L_j, L_p, L_i, L_m

      real(rk) :: K
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type type_population

   ! Standard variable ("total mass") used for mass conservation checking
   type (type_bulk_standard_variable),parameter :: total_energy = type_bulk_standard_variable(name='total_energy',units='J',aggregate_variable=.true.,conserved=.true.)
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the module
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
   class (type_population), intent(inout),target :: self
   integer,                                intent(in )          :: configunit
!
! !LOCAL VARIABLES:
   integer :: iclass
   real(rk) :: delta_logV
   real(rk) :: V_min, V_max
   character(len=8) :: strindex
   real(rk) :: NV_ini, NE_ini
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read parameters
   ! All rate coefficients are converted from d-1 to s-1 ensure source terms are directly compatible with FABM.
   ! Default values taken from Table S5
   call self%get_parameter(self%nclass,'nclass', '',     'number of structural volume classes', default=100)
   call self%get_parameter(self%p_Am, 'p_Am',   'J cm-2 d-1', 'maximum specific assimilation rate')
   call self%get_parameter(self%p_T,  'p_T',    'J cm-2 d-1', 'maximum surface-area-specific maintenance rate')
   call self%get_parameter(self%p_M,  'p_M',    'J cm-3 d-1', 'maximum specific assimilation rate')
   call self%get_parameter(self%E_Hb, 'E_Hb',   'J',          'maturity at birth')
   call self%get_parameter(self%E_Hj, 'E_Hj',   'J',          'maturity at metamorphosis')
   call self%get_parameter(self%E_Hp, 'E_Hp',   'J',          'maturity at puberty')
   call self%get_parameter(self%v,    'v',      'cm d-1',     'energy conductance')
   call self%get_parameter(self%k_J,  'k_J',    'd-1',        'maturity maintenance rate')
   call self%get_parameter(self%kap,  'kap',    '-',          'investment in soma (remainder goes to reproduction)')
   call self%get_parameter(self%kap_R,'kap_R',  '-',          'reproduction efficiency')
   call self%get_parameter(self%E_G,  'E_G',    'J cm-3',     'energy density of structure')
   call self%get_parameter(self%h_a,  'h_a',    'd-2',        'ageing acceleration')
   call self%get_parameter(self%s_G,  's_G',    '-',          'stress coefficient')
   call self%get_parameter(self%K,    'K',      'J',          'food half saturation')

   ! for now: get crucial implied properties as parameter rather than self computing them
   call self%get_parameter(self%E_0,  'E_0', 'J',  'initial energy content of an egg')
   call self%get_parameter(self%L_b,  'L_b', 'cm', 'structural length at birth')
   call self%get_parameter(self%L_j,  'L_j', 'cm', 'structural length at metamorphosis')
   self%s_M = self%L_j/self%L_b

   self%L_m = self%kap*self%p_Am/self%p_M
   self%L_T = self%p_T/self%p_M
   self%L_i = (self%L_m - self%L_T)*self%s_M
   self%E_m = self%p_Am/self%v

   ! Determine size classes (log-spaced between size at birth and infinite size)
   allocate(self%logV_center(self%nclass))
   allocate(self%V_center(self%nclass))
   allocate(self%L_center(self%nclass))
   allocate(self%delta_V(self%nclass))
   V_min = 0.1_rk*self%L_b**3
   V_max = self%L_i**3
   delta_logV = (log(V_max)-log(V_min))/(self%nclass-1)
   do iclass=1,self%nclass
      self%logV_center(iclass) = log(V_min)+delta_logV*(iclass-1)
   end do
   self%V_center = exp(self%logV_center)                                         ! Mass of each size class
   self%L_center = self%V_center**(1._rk/3._rk)
   self%delta_V = exp(self%logV_center+delta_logV) - self%V_center               ! Mass difference between consecutive size classes

   allocate(self%id_NV(self%nclass))
   allocate(self%id_NE(self%nclass))
   allocate(self%id_NE_H(self%nclass))
   allocate(self%id_NQ(self%nclass))
   allocate(self%id_Nh(self%nclass))
   do iclass=1,self%nclass
      write (strindex, '(i0)') iclass
      if (iclass == 1) then
         NV_ini = self%V_center(iclass)
         NE_ini = self%E_0
      else
         NV_ini = 0
         NE_ini = 0
      end if
      call self%register_state_variable(self%id_NV(iclass), 'NV_'//trim(strindex), 'cm3', 'structural volume in bin '//trim(strindex), initial_value=NV_ini)
      call self%register_state_variable(self%id_NE(iclass), 'NE_'//trim(strindex), 'J', 'reserve in bin '//trim(strindex), initial_value=NE_ini)
      call self%register_state_variable(self%id_NE_H(iclass), 'NE_H_'//trim(strindex), 'J', 'maturity in bin '//trim(strindex))
      call self%register_state_variable(self%id_NQ(iclass), 'NQ_'//trim(strindex), '-', 'damage inducing compounds in bin '//trim(strindex))
      call self%register_state_variable(self%id_Nh(iclass), 'Nh_'//trim(strindex), 'd-1 cm3', 'hazard in bin '//trim(strindex))
      call self%set_variable_property(self%id_NV(iclass),'V',self%V_center(iclass))
      call self%add_to_aggregate_variable(total_energy, self%id_NV(iclass), scale_factor=self%E_G)
      call self%add_to_aggregate_variable(total_energy, self%id_NE(iclass))
   end do
   call self%register_state_variable(self%id_waste, 'waste', 'J', 'waste target')
   call self%add_to_aggregate_variable(total_energy, self%id_waste)

   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)

   call self%register_diagnostic_variable(self%id_R, 'R', '# d-1', 'reproduction rate')
   call self%register_diagnostic_variable(self%id_N_e, 'N_e', '#', 'number of non-feeding individuals (eggs, embryos)')
   call self%register_diagnostic_variable(self%id_N_j, 'N_j', '#', 'number of juveniles (feeding but not reproducing)')
   call self%register_diagnostic_variable(self%id_N_a, 'N_a', '#', 'number of adults (feeding and reproducing)')

   self%dt = 86400

   end subroutine initialize
!EOC

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_population),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: v_E_G_plus_p_T_per_kap
      real(rk) :: p_M_per_kap
      real(rk) :: p_T_per_kap
      real(rk) :: E_G_per_kap
      real(rk) :: one_minus_kap
      real(rk) :: L_m3
      real(rk) :: waste
      
      integer :: iclass
      real(rk) :: f
      real(rk) :: L, L2, L3, s, p_A, p_C, p_R, p_R_sum, R_p, dV
      real(rk) :: N_e, N_j, N_a
      real(rk), dimension(self%nclass) :: NV, NE, NE_H, NQ, Nh
      real(rk), dimension(self%nclass) :: N, hazard
      real(rk), dimension(self%nclass) :: dE, dL, dE_H, dQ, dH
      real(rk), dimension(0:self%nclass) :: nflux, Eflux, E_Hflux, Qflux, hflux, E, E_H, Q, H
      real(rk), parameter :: delta_t = 12._rk/86400

      v_E_G_plus_p_T_per_kap = (self%v*self%E_G + self%p_T)/self%kap
      p_M_per_kap = self%p_M/self%kap
      p_T_per_kap = self%p_T/self%kap
      E_G_per_kap = self%E_G/self%kap
      one_minus_kap = 1.0_rk - self%kap
      L_m3 = self%L_m**3
      
      _HORIZONTAL_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_waste, waste)
         N_e = 0
         N_j = 0
         N_a = 0
         do iclass=1,self%nclass
            _GET_HORIZONTAL_(self%id_NV(iclass), NV(iclass))
            _GET_HORIZONTAL_(self%id_NE(iclass), NE(iclass))
            _GET_HORIZONTAL_(self%id_NE_H(iclass), NE_H(iclass))
            _GET_HORIZONTAL_(self%id_NQ(iclass), NQ(iclass))
            _GET_HORIZONTAL_(self%id_NH(iclass), Nh(iclass))
      
            N(iclass) = max(0.0_rk, NV(iclass))/self%V_center(iclass)
            if (N(iclass) > 0) then
               E(iclass) = max(0.0_rk, NE(iclass)/N(iclass))
               E_H(iclass) = max(0.0_rk, NE_H(iclass)/N(iclass))
               Q(iclass) = max(0.0_rk, NQ(iclass)/N(iclass))
               H(iclass) = max(0.0_rk, Nh(iclass)/N(iclass))
               if (E_H(iclass) > self%E_Hp) then
                  N_a = N_a + N(iclass)
               elseif (E_H(iclass) > self%E_Hb) then
                  N_j = N_j + N(iclass)
               else
                  N_e = N_e + N(iclass)
               end if
            else
               E(iclass) = 0
               E_H(iclass) = 0
               Q(iclass) = 0
               H(iclass) = 0
            end if
         end do
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_N_e, N_e)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_N_j, N_j)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_N_a, N_a)
         hazard = h/self%V_center

         ! Individual physiology (per size class)
         p_R_sum = 0
         dE_H = 0
         f = max(0._rk, waste)/(max(0._rk, waste) + self%K)
         do iclass=1,self%nclass
            L = self%L_center(iclass)
            L3 = self%V_center(iclass)
            L2 = L*L
            s = max(1.0_rk, min(self%s_M, L/self%L_b))

            ! Energy fluxes in J/d per individual
            if (E_H(iclass) < self%E_Hb) then
               ! Non-feeding embryo
               p_A = 0
            else
               ! Feeding juvenile or adult
               p_A = self%p_Am*L2*f*s
            end if
            p_C = E(iclass)*(v_E_G_plus_P_T_per_kap*s*L2 + p_M_per_kap*L3)/(E(iclass) + E_G_per_kap*L3)
            p_R = one_minus_kap*p_C - self%k_J*E_H(iclass)

            ! Change in reserve (J), structural length (cm), maturity (J), reproduction buffer (J)
            dE(iclass) = p_A - p_C
            dL(iclass) = (E(iclass)*self%v*s-(p_M_per_kap*L+p_T_per_kap*s)*L3)/3/(E(iclass)+E_G_per_kap*L3)

            ! Verify consistence between expressions for j_C and r
            !write (*,*) p_C, (self%v/L - dL(iclass)/L*3)*E(iclass), p_C - (self%v/L - dL(iclass)/L*3)*E(iclass)
            !write (*,*) (self%kap*p_C - self%p_M*L3 - self%p_T*L2)/self%E_G, dL(iclass)*3*L2, (self%kap*p_C - self%p_M*L3 - self%p_T*L2)/self%E_G - dL(iclass)*3*L2

            if (E_H(iclass) < self%E_Hp) then
                ! Embryo or juvenile allocating to maturation
                dE_H(iclass) = p_R
            else
                ! Mature adult allocating to reproduction
                p_R_sum = p_R_sum + N(iclass)*p_R
            end if

            ! Damage-inducing compounds, damage, survival (0-1) - p 216
            dQ(iclass) = (Q(iclass)/L_m3*self%s_G + self%h_a)*max(0., p_C)/self%E_m
            dH(iclass) = Q(iclass)

            _SET_BOTTOM_ODE_(self%id_waste, N(iclass)*(-p_A + self%p_M*L3 + self%p_T*L2 + one_minus_kap*p_C))
         end do

         ! Compute number of individuals moving from each size class to the next (units: # s-1)
         nflux = 0
         do iclass=1,self%nclass
            dV = 3*dL(iclass)*self%L_center(iclass)**2
            if (dV >= 0) then
               ! Growing: positive flux (towards larger size) over right boundary
               nflux(iclass) = nflux(iclass) + N(iclass)*dV/self%delta_V(iclass)
            else
               ! Growing: negative flux (towards smaller size) over left boundary
               nflux(iclass-1) = nflux(iclass-1) + N(iclass)*dV/self%delta_V(iclass-1)
            end if
         end do

         R_p = self%kap_R*p_R_sum/self%E_0/2  ! Divide by 2 assuming sexual reproduction (only female investment counts)
         _SET_BOTTOM_ODE_(self%id_waste, -R_p*self%E_0)

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R, R_p)

         ! Use recruitment as number of incoming individuals for the first size class.
         nflux(0) = R_p
         E(0) = self%E_0 - self%V_center(1)*self%E_G
         E_H(0) = 0
         Q(0) = 0
         h(0) = 0

         nflux(self%nclass) = 0
         do iclass=0,self%nclass
            if (nflux(iclass) >= 0) then
               ! Positive flux of individuals (towards larger size)
               Eflux(iclass) = nflux(iclass)*E(iclass)
               E_Hflux(iclass) = nflux(iclass)*E_H(iclass)
               Qflux(iclass) = nflux(iclass)*Q(iclass)
               hflux(iclass) = nflux(iclass)*h(iclass)
            else
               ! Negative flux of individuals (towards smaller size)
               Eflux(iclass) = nflux(iclass)*E(iclass+1)
               E_Hflux(iclass) = nflux(iclass)*E_H(iclass+1)
               Qflux(iclass) = nflux(iclass)*Q(iclass+1)
               hflux(iclass) = nflux(iclass)*h(iclass+1)
            end if
         end do
         
         ! Transfer size-class-specific source terms and diagnostics to FABM
         do iclass=1,self%nclass
            ! Apply specific mortality (s-1) to size-class-specific abundances and apply upwind advection - this is a time-explicit version of Eq G.1 of Hartvig et al.
            _SET_BOTTOM_ODE_(self%id_NV(iclass),  -hazard(iclass)*NV(iclass)   + (nflux (iclass-1) - nflux  (iclass))*self%V_center(iclass))
            _SET_BOTTOM_ODE_(self%id_NE(iclass),  -hazard(iclass)*NE(iclass)   + Eflux  (iclass-1) - Eflux  (iclass) + dE  (iclass)*N(iclass))
            _SET_BOTTOM_ODE_(self%id_NE_H(iclass),-hazard(iclass)*NE_H(iclass) + E_Hflux(iclass-1) - E_Hflux(iclass) + dE_H(iclass)*N(iclass))
            _SET_BOTTOM_ODE_(self%id_NQ(iclass),  -hazard(iclass)*NQ(iclass)   + Qflux  (iclass-1) - Qflux  (iclass) + dQ  (iclass)*N(iclass))
            _SET_BOTTOM_ODE_(self%id_Nh(iclass),  -hazard(iclass)*Nh(iclass)   + hflux  (iclass-1) - Hflux  (iclass) + dH  (iclass)*N(iclass))
         end do
         _SET_BOTTOM_ODE_(self%id_waste,sum(hazard*(NV*self%E_G + NE)))
      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

!-----------------------------------------------------------------------

end module deb_population

!-----------------------------------------------------------------------
! Copyright Jorn Bruggeman/PML 2015-2017
!-----------------------------------------------------------------------
