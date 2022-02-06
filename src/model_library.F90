module deb_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use deb_population
   ! Add new DEB models here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: deb_model_factory

contains

   subroutine create(self,name,model)

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('population'); allocate(type_population::model)
         ! Add new DEB models here
      end select

   end subroutine

end module
