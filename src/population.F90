#include "fabm_driver.h"

! DEB structured population tailored to Artemia urmiana
!
! This model structures the population along a structural volume axis,
! which can be partitioned over an arbitrary number of bins (size classes).
! The biomass per size class is represented by a depth-integrated field.
!
! In addition, each size class has associated reserve (E), maturity plus,
! reproduction buffer (E_HR), damage-inducing compounds (Q) and damage (H, equivalent
! to structure-weighted hazard). These properties are represented by the totals for a
! single size class. That is, they are summed over all individuals in the size class.
! For instance, reserve in that size class is expressed as NE, which equals the
! number of individuals in that size class times their mean reserve content in Joule.
! This representation ensures mass conservation.
!
! Temperature tolerance is described by an Arrhenius relationship.
!
! Energy allocated to reproduction is placed into a buffer that is converted into offspring,
! at some prescribed, temperature-dependent rate.
!
! For the purpose of describing Artemia, the maintenance rate has been made dependent on salinity.
! Currently a DEBtox-like functional relationship is used for this.

module deb_population

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_population
      ! Variable identifiers
      type (type_bottom_state_variable_id), allocatable :: id_NV(:)
      type (type_bottom_state_variable_id), allocatable :: id_NE(:)
      type (type_bottom_state_variable_id), allocatable :: id_NE_HR(:)
      type (type_bottom_state_variable_id), allocatable :: id_NQ(:)
      type (type_bottom_state_variable_id), allocatable :: id_NH(:)

      type (type_bottom_state_variable_id) :: id_waste, id_food

      type (type_bottom_diagnostic_variable_id) :: id_R, id_E_R_sum
      type (type_bottom_diagnostic_variable_id) :: id_N_b, id_N_j, id_N_p, id_N_a
      type (type_bottom_diagnostic_variable_id) :: id_Et_b, id_Et_j, id_Et_p, id_Et_a

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

      real(rk) :: T_A
      real(rk) :: T_ref

      real(rk) :: K

      real(rk) :: R_bg

      real(rk) :: salt_threshold
      real(rk) :: p_T_at_s300

      real(rk) :: reproduction_frequency
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type type_population

   ! Standard variable ("total mass") used for mass conservation checking
   type (type_bulk_standard_variable),parameter :: total_energy = type_bulk_standard_variable(name='total_energy',units='J m-3',aggregate_variable=.true.,conserved=.true.)

contains

   subroutine initialize(self,configunit)
   class (type_population), intent(inout),target :: self
   integer,                                intent(in )          :: configunit

   integer :: iclass
   real(rk) :: delta_logV
   real(rk) :: V_min, V_max
   character(len=8) :: strindex
   real(rk) :: N_ini

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

   ! Temperature dependence
   call self%get_parameter(self%T_A,  'T_A',    'K',          'Arrhenius temperature (activation energy divided by universal gas constant)')
   call self%get_parameter(self%T_ref,'T_ref',  'degrees_C',  'Reference temperature', default=20._rk)

   ! Salinity dependence
   call self%get_parameter(self%salt_threshold,   'salt_threshold',   '-', 'salinity above which surface-area-specific maintenance rate starts to increase', default=0._rk)
   call self%get_parameter(self%p_T_at_s300, 'p_T_at_s300', '-', 'surface-area-specific maintenance rate at 300 PSU', default=self%p_T)

   ! Half-saturation for food
   call self%get_parameter(self%K,    'K',      'J',          'food half saturation')

   call self%get_parameter(self%reproduction_frequency, 'reproduction_frequency', 'd-1', 'reproduction frequency', default=0.1_rk)
   call self%get_parameter(self%R_bg,         'R_bg',         '# d-1', 'emergence from non-tracked resting eggs (cysts)', default=0.0_rk)

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
   V_min = 0.01_rk*self%L_b**3
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
   allocate(self%id_NE_HR(self%nclass))
   allocate(self%id_NQ(self%nclass))
   allocate(self%id_Nh(self%nclass))
   do iclass=1,self%nclass
      write (strindex, '(i0)') iclass
      if (iclass == 1) then
         N_ini = 1._rk
      else
         N_ini = 0._rk
      end if
      call self%register_state_variable(self%id_NV(iclass), 'NV_'//trim(strindex), 'cm3 m-2', 'structural volume in bin '//trim(strindex), initial_value=N_ini*self%V_center(iclass))
      call self%register_state_variable(self%id_NE(iclass), 'NE_'//trim(strindex), 'J m-2', 'reserve in bin '//trim(strindex), initial_value=N_ini*(self%E_0 - self%V_center(iclass)*self%E_G))
      call self%register_state_variable(self%id_NE_HR(iclass), 'NE_HR_'//trim(strindex), 'J m-2', 'maturity + reproduction buffer in bin '//trim(strindex))
      call self%register_state_variable(self%id_NQ(iclass), 'NQ_'//trim(strindex), 'm-2', 'damage-inducing compounds in bin '//trim(strindex))
      call self%register_state_variable(self%id_Nh(iclass), 'NH_'//trim(strindex), 'd-1 cm3 m-2', 'structure-weighted hazard in bin '//trim(strindex))
      call self%set_variable_property(self%id_NV(iclass),'V',self%V_center(iclass))
      call self%add_to_aggregate_variable(total_energy, self%id_NV(iclass), scale_factor=self%E_G)
      call self%add_to_aggregate_variable(total_energy, self%id_NE(iclass))
   end do
   call self%register_state_variable(self%id_waste, 'waste', 'J m-2', 'waste target')
   call self%register_state_dependency(self%id_food, 'food', 'J m-2', 'food source')
   call self%add_to_aggregate_variable(total_energy, self%id_waste)

   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)

   call self%register_diagnostic_variable(self%id_R, 'R', '# m-2 d-1', 'reproduction rate')
   call self%register_diagnostic_variable(self%id_N_b, 'N_b', '# m-2', 'number of eggs')
   call self%register_diagnostic_variable(self%id_N_j, 'N_j', '# m-2', 'number of juveniles pre-metamorphosis')
   call self%register_diagnostic_variable(self%id_N_p, 'N_p', '# m-2', 'number of juveniles post-metamorphosis')
   call self%register_diagnostic_variable(self%id_N_a, 'N_a', '# m-2', 'number of adults')
   call self%register_diagnostic_variable(self%id_Et_b, 'Et_b', 'J m-2', 'energy in eggs')
   call self%register_diagnostic_variable(self%id_Et_j, 'Et_j', 'J m-2', 'energy in juveniles pre-metamorphosis')
   call self%register_diagnostic_variable(self%id_Et_p, 'Et_p', 'J m-2', 'energy in juveniles post-metamorphosis')
   call self%register_diagnostic_variable(self%id_Et_a, 'Et_a', 'J m-2', 'energy in adults')
   call self%register_diagnostic_variable(self%id_E_R_sum, 'E_R_sum', 'J m-2', 'total reproduction buffer')
   call self%add_to_aggregate_variable(total_energy, self%id_E_R_sum)

   self%dt = 86400

   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_population),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: v_E_G_plus_p_T_per_kap
      real(rk) :: p_M_per_kap
      real(rk) :: p_T_per_kap
      real(rk) :: E_G_per_kap
      real(rk) :: one_minus_kap
      real(rk) :: food, temp, salt, T_K
      real(rk) :: c_T
      real(rk) :: v, k_J, p_Am, p_M, p_T, h_a
      real(rk) :: s_G_per_L_m3_E_m, h_a_per_E_m
      real(rk) :: invdenom

      integer :: iclass
      real(rk) :: f
      real(rk) :: L, L2, L3, s, p_A, p_C, p_R_sum, R_p, dV, E_R_sum
      real(rk) :: N_b, N_j, N_p, N_a
      real(rk) :: Et_b, Et_j, Et_p, Et_a
      real(rk), dimension(self%nclass) :: NV, NE, NE_HR, NQ, NH
      real(rk), dimension(self%nclass) :: N, hazard, E_H, E_R
      real(rk), dimension(self%nclass) :: dE, dL, dE_HR, dQ, dH
      real(rk), dimension(0:self%nclass) :: nflux, Eflux, E_HRflux, Qflux, hflux, E, E_HR, Q, H
      real(rk), parameter :: delta_t = 12._rk/86400
      real(rk), parameter :: Kelvin = 273.15_rk
      real(rk), parameter :: N_eps = 1e-12_rk
      real(rk), parameter :: onethird = 1._rk / 3._rk

      _BOTTOM_LOOP_BEGIN_

         _GET_HORIZONTAL_(self%id_food, food)

         _GET_(self%id_temp, temp)
         _GET_(self%id_salt, salt)
         salt = max(0._rk, salt)

         ! Apply temperature dependence of rates
         T_K = temp + Kelvin
         c_T = exp(self%T_A / (self%T_ref + Kelvin) - self%T_A / T_K)
         v = self%v * c_T
         k_J = self%k_J * c_T
         p_Am = self%p_Am * c_T
         p_M = self%p_M * c_T
         p_T = (self%p_T + max(salt - self%salt_threshold, 0.) * (self%p_T_at_s300 - self%p_T) / (300._rk - self%salt_threshold)) * c_T
         h_a = self%h_a * c_T * c_T

         E_G_per_kap = self%E_G / self%kap
         p_M_per_kap = p_M / self%kap
         p_T_per_kap = p_T / self%kap
         v_E_G_plus_p_T_per_kap = (v * self%E_G + p_T) / self%kap
         s_G_per_L_m3_E_m = self%s_G / self%L_m**3 / self%E_m
         h_a_per_E_m = h_a / self%E_m
         one_minus_kap = 1.0_rk - self%kap

         do iclass=1,self%nclass
            _GET_HORIZONTAL_(self%id_NV(iclass), NV(iclass))         ! total structural volume in bin (cm3)
            _GET_HORIZONTAL_(self%id_NE(iclass), NE(iclass))         ! total reserve in bin (J)
            _GET_HORIZONTAL_(self%id_NE_HR(iclass), NE_HR(iclass))   ! total energy allocated to maturity and reproduction in bin (J)
            _GET_HORIZONTAL_(self%id_NQ(iclass), NQ(iclass))         ! damage (d-2) times total structural biomass (cm3) in bin
            _GET_HORIZONTAL_(self%id_NH(iclass), Nh(iclass))         ! hazard rate (d-1) times total structural biomass (cm3) in bin

            N(iclass) = max(0.0_rk, NV(iclass)) / self%V_center(iclass)  ! compute number of individuals in bin
            if (N(iclass) > N_eps) then
               ! Positive number of individuals in this class
               E(iclass) = max(0.0_rk, NE(iclass) / N(iclass))         ! reserve per individual (J)
               E_HR(iclass) = max(0.0_rk, NE_HR(iclass) / N(iclass))   ! energy allocated to maturity and reproduction per individual (J)
               E_H(iclass) = min(E_HR(iclass), self%E_Hp)              ! energy allocated to maturity per individual (J)
               E_R(iclass) = max(0.0_rk, E_HR(iclass) - self%E_Hp)     ! energy allocated to reproduction per individual (J)
               Q(iclass) = max(0.0_rk, NQ(iclass) / N(iclass))         ! damage (d-2) times structural biomass (cm3) per individual
               H(iclass) = max(0.0_rk, NH(iclass) / N(iclass))         ! hazard rate (d-1) times total structural biomass (cm3) per individual
            else
               ! Non-positive number of individuals in this class
               E(iclass) = 0
               E_HR(iclass) = 0
               E_H(iclass) = 0
               E_R(iclass) = 0
               Q(iclass) = 0
               H(iclass) = 0
            end if
         end do
         E_R_sum = sum(E_R(1:self%nclass) * N(1:self%nclass))
         _SET_BOTTOM_DIAGNOSTIC_(self%id_E_R_sum, E_R_sum)   ! total energy allocated to reproduction (for mass conservation checking)

         ! Count individuals per life stage
         N_b = 0
         N_j = 0
         N_p = 0
         N_a = 0
         Et_b = 0
         Et_j = 0
         Et_p = 0
         Et_a = 0
         do iclass=1,self%nclass
            if (E_H(iclass) == self%E_Hp) then
               ! Adult
               N_a = N_a + N(iclass)
               Et_a = Et_a + NV(iclass) * self%E_G + NE(iclass)
            elseif (E_H(iclass) > self%E_Hj) then
               ! Juvenile post metamorphosis
               N_p = N_p + N(iclass)
               Et_p = Et_p + NV(iclass) * self%E_G + NE(iclass)
            elseif (E_H(iclass) > self%E_Hb) then
               ! Juvenile pre metamorphosis
               N_j = N_j + N(iclass)
               Et_j = Et_j + NV(iclass) * self%E_G + NE(iclass)
            else
               ! Egg/embryo
               N_b = N_b + N(iclass)
               Et_b = Et_b + NV(iclass) * self%E_G + NE(iclass)
            end if
         end do
         _SET_BOTTOM_DIAGNOSTIC_(self%id_N_b, N_b)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_N_j, N_j)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_N_p, N_p)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_N_a, N_a)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_Et_b, Et_b)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_Et_j, Et_j)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_Et_p, Et_p)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_Et_a, Et_a)

         hazard = H(1:self%nclass) / self%V_center

         ! Individual physiology (per size class)
         dE_HR = 0
         f = max(0._rk, food)/(max(0._rk, food) + self%K)
         !f = 1
         do iclass=1,self%nclass
            L = self%L_center(iclass)
            L3 = self%V_center(iclass)
            L2 = L*L
            s = max(1.0_rk, min(self%s_M, L/self%L_b))

            ! Energy fluxes in J/d per individual
            if (E_H(iclass) < self%E_Hb) then
               ! Non-feeding embryo
               p_A = 0._rk
            else
               ! Feeding juvenile or adult
               p_A = p_Am * f * L2 * s
            end if

            ! Catabolic flux p_C (J/d)
            invdenom = 1.0_rk / (E(iclass) + E_G_per_kap * L3)
            p_C = E(iclass) * (v_E_G_plus_P_T_per_kap * s * L2 + p_M_per_kap * L3) * invdenom

            ! Change in reserve E (J) and structural length L (cm)
            dE(iclass) = p_A - p_C
            dL(iclass) = (E(iclass) * v * s - (p_M_per_kap * L + p_T_per_kap * s) * L3) * onethird * invdenom

            ! Allocation to maturity/reproduction (J/d) (maturity maintenance already subtracted)
            dE_HR(iclass) = one_minus_kap * p_C - k_J * E_H(iclass)

            ! Verify consistence between expressions for j_C and r
            !write (*,*) p_C, (v/L - dL(iclass)/L*3)*E(iclass), p_C - (v/L - dL(iclass)/L*3)*E(iclass)
            !write (*,*) (self%kap*p_C - p_M*L3 - p_T*L2)/self%E_G, dL(iclass)*3*L2, (self%kap*p_C - p_M*L3 - p_T*L2)/self%E_G - dL(iclass)*3*L2

            ! Change in damage-inducing compounds and structure-weighted hazard - p 216 in DEB 2010
            dQ(iclass) = (Q(iclass) * s_G_per_L_m3_E_m + h_a_per_E_m) * max(0., p_C)
            dH(iclass) = Q(iclass)

            ! Take assimilated energy away from food source.
            _ADD_BOTTOM_SOURCE_(self%id_food, -N(iclass) * p_A)

            ! Add energy fluxes towards maintenance and maturity
            _ADD_BOTTOM_SOURCE_(self%id_waste, N(iclass) * (p_M * L3 + p_T * L2 * s + k_J * E_H(iclass)))
            if (E_H(iclass) < self%E_Hp) _ADD_BOTTOM_SOURCE_(self%id_waste,  N(iclass) * dE_HR(iclass))
         end do

         ! Compute number of individuals moving from each size class to the next (units: # d-1)
         nflux = 0
         do iclass=1,self%nclass
            ! Compute change in structural volume V from change in structural length L (V=L^3)
            ! NB the change in V is not well-defined [always zero] at V=0 (conception) - unlike the change in L!
            dV = 3 * dL(iclass) * self%L_center(iclass)**2
            if (dV >= 0) then
               ! Growing: positive flux of individuals (towards larger size) over right boundary
               nflux(iclass) = nflux(iclass) + N(iclass) * dV
            else
               ! Shrinking: negative flux of individuals (towards smaller size) over left boundary
               nflux(iclass-1) = nflux(iclass-1) + N(iclass) * dV
            end if
         end do

         ! Divide structural volume flux between bins by distance between bin centers to arive at fluxes of individuals.
         nflux(1:self%nclass - 1) = nflux(1:self%nclass - 1) / self%delta_V(1:self%nclass - 1)

         ! Energy is first placed in reproduction buffer, then converted into offspring at a particular rate.
         p_R_sum = c_T * self%reproduction_frequency * E_R_sum
         do iclass=1,self%nclass
            dE_HR(iclass) = dE_HR(iclass) - c_T * self%reproduction_frequency * E_R(iclass)
         end do

         ! Compute reproduction (# d-1) from population-integrated energy flux allocated to reproduction.
         ! NB we divide by 2 assuming sexual reproduction (only female investment leads to eggs)
         R_p = 0.5_rk * self%kap_R * p_R_sum / self%E_0
         _SET_BOTTOM_DIAGNOSTIC_(self%id_R, R_p)

         ! Emergence from non-tracked resting cysts [parameter is given at reference temperature]
         R_p = R_p + self%R_bg * c_T

         ! Energy lost during reproduction goes to waste
         _ADD_BOTTOM_SOURCE_(self%id_waste, p_R_sum - R_p * self%E_0)

         ! Use recruitment as number of incoming individuals for the first size class.
         ! Reserve per individual is the reserve per egg, minus the energy contained in structure at the lowest tracked size.
         nflux(0) = R_p
         E(0) = self%E_0 - self%V_center(1) * self%E_G
         E_HR(0) = 0
         Q(0) = 0
         h(0) = 0

         ! Compute fluxes of reserve, maturity, damange-inducing compounds, structure-weighted hazard between size classes.
         nflux(self%nclass) = 0
         do iclass=0,self%nclass
            if (nflux(iclass) >= 0) then
               ! Positive flux of individuals (towards larger size)
               Eflux   (iclass) = nflux(iclass) * E   (iclass)
               E_HRflux(iclass) = nflux(iclass) * E_HR(iclass)
               Qflux   (iclass) = nflux(iclass) * Q   (iclass)
               hflux   (iclass) = nflux(iclass) * h   (iclass)
            else
               ! Negative flux of individuals (towards smaller size)
               Eflux   (iclass) = nflux(iclass) * E   (iclass+1)
               E_HRflux(iclass) = nflux(iclass) * E_HR(iclass+1)
               Qflux   (iclass) = nflux(iclass) * Q   (iclass+1)
               hflux   (iclass) = nflux(iclass) * h   (iclass+1)
            end if
         end do

         ! Transfer size-class-specific source terms to FABM
         do iclass=1,self%nclass
            ! Apply specific mortality (s-1) to size-class-specific abundances and apply upwind advection - this is a time-explicit version of Eq G.1 of Hartvig et al.
            _ADD_BOTTOM_SOURCE_(self%id_NV(iclass),   -hazard(iclass)*NV   (iclass) + (nflux  (iclass-1) - nflux   (iclass))*self%V_center(iclass))
            _ADD_BOTTOM_SOURCE_(self%id_NE(iclass),   -hazard(iclass)*NE   (iclass) + Eflux   (iclass-1) - Eflux   (iclass) + dE   (iclass)*N(iclass))
            _ADD_BOTTOM_SOURCE_(self%id_NE_HR(iclass),-hazard(iclass)*NE_HR(iclass) + E_HRflux(iclass-1) - E_HRflux(iclass) + dE_HR(iclass)*N(iclass))
            _ADD_BOTTOM_SOURCE_(self%id_NQ(iclass),   -hazard(iclass)*NQ   (iclass) + Qflux   (iclass-1) - Qflux   (iclass) + dQ   (iclass)*N(iclass))
            _ADD_BOTTOM_SOURCE_(self%id_Nh(iclass),   -hazard(iclass)*NH   (iclass) + hflux   (iclass-1) - Hflux   (iclass) + dH   (iclass)*N(iclass))
         end do

         ! Add dead structure and reserve to waste pool.
         _ADD_BOTTOM_SOURCE_(self%id_waste, sum(hazard(1:self%nclass) * (NV(1:self%nclass) * self%E_G + NE(1:self%nclass) + E_R(1:self%nclass) * N(1:self%nclass))))

      _BOTTOM_LOOP_END_

   end subroutine do_bottom

!-----------------------------------------------------------------------

end module deb_population

!-----------------------------------------------------------------------
! Copyright Jorn Bruggeman/PML 2015-2017
!-----------------------------------------------------------------------
