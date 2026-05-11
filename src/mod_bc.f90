!> Boundary condition management and evaluation.
!!
!! This module handles the mapping of case parameters to mesh patches
!! and provides routines for evaluating boundary values for velocity,
!! pressure, and species.
module mod_bc
   use mod_kinds, only : rk, zero, name_len, fatal_error, lowercase
   use mod_input, only : case_params_t, max_species
   use mod_mesh_types, only : mesh_t
   implicit none

   private

   integer, parameter, public :: bc_unknown   = 0
   integer, parameter, public :: bc_wall      = 1
   integer, parameter, public :: bc_symmetry  = 2
   integer, parameter, public :: bc_periodic  = 3
   integer, parameter, public :: bc_dirichlet = 4
   integer, parameter, public :: bc_neumann   = 5

   !> Individual boundary patch settings.
   type, public :: bc_patch_t
      integer :: patch_id = 0                    !< Patch ID from mesh
      character(len=name_len) :: name = ""       !< Patch name
      character(len=name_len) :: type_name = ""  !< Boundary type name (e.g., "wall")
      integer :: type_id = bc_unknown            !< Internal type ID
      real(rk) :: velocity(3) = zero             !< Specified velocity vector
      real(rk) :: pressure = zero                !< Specified pressure value
      real(rk) :: dpdn = zero                    !< Specified pressure gradient

      !> Species boundary conditions
      integer :: species_type_id = bc_unknown    !< BC type for species
      real(rk) :: species_Y(max_species) = zero  !< Specified mass fractions
   end type bc_patch_t

   !> Complete set of boundary conditions for all patches.
   type :: bc_set_t
      integer :: npatches = 0                        !< Number of patches
      type(bc_patch_t), allocatable :: patches(:)    !< Array of patch BCs
   end type bc_set_t

   public :: build_bc_set, finalize_bc_set, bc_set_t
   public :: patch_type_for_face, boundary_velocity, face_effective_neighbor
   public :: boundary_species
   public :: is_periodic_face

contains

   !> Map case parameters to mesh patches to create the BC set.
   !!
   !! @param mesh Mesh data structure.
   !! @param params Case configuration parameters.
   !! @param bc Boundary condition set to build.
   subroutine build_bc_set(mesh, params, bc)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(bc_set_t), intent(inout) :: bc

      integer :: p, q
      logical :: found

      call finalize_bc_set(bc)
      bc%npatches = mesh%npatches
      allocate(bc%patches(mesh%npatches))

      do p = 1, mesh%npatches
         bc%patches(p)%patch_id = p
         bc%patches(p)%name = mesh%patches(p)%name
         bc%patches(p)%type_name = "wall"
         bc%patches(p)%type_id = bc_wall
         found = .false.

         do q = 1, params%n_patches
            if (trim(params%patch_name(q)) == trim(mesh%patches(p)%name)) then
               bc%patches(p)%type_name = trim(lowercase(params%patch_type(q)))
               bc%patches(p)%type_id = parse_bc_type(params%patch_type(q))
               bc%patches(p)%velocity = [params%patch_u(q), params%patch_v(q), params%patch_w(q)]
               bc%patches(p)%pressure = params%patch_p(q)
               bc%patches(p)%dpdn = params%patch_dpdn(q)

               if (trim(params%patch_species_type(q)) /= "") then
                  bc%patches(p)%species_type_id = parse_bc_type(params%patch_species_type(q))
               else
                  bc%patches(p)%species_type_id = parse_bc_type(params%patch_type(q))
               end if
               bc%patches(p)%species_Y(:) = params%patch_Y(:, q)

               found = .true.
               exit
            end if
         end do

         if (.not. found) then
            call fatal_error('bc', 'mesh patch '//trim(mesh%patches(p)%name)//' missing from case.nml')
         end if

         if (bc%patches(p)%type_id == bc_periodic) then
            call ensure_periodic_links(mesh, p)
         end if
      end do
   end subroutine build_bc_set


   !> Deallocate the boundary condition set.
   !!
   !! @param bc Boundary condition set to finalize.
   subroutine finalize_bc_set(bc)
      type(bc_set_t), intent(inout) :: bc

      if (allocated(bc%patches)) deallocate(bc%patches)
      bc%npatches = 0
   end subroutine finalize_bc_set


   integer function parse_bc_type(text) result(type_id)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: key

      type_id = bc_unknown

      key = trim(lowercase(text))
      select case (trim(key))
      case ('wall')
         type_id = bc_wall
      case ('symmetry', 'symmetric')
         type_id = bc_symmetry
      case ('periodic')
         type_id = bc_periodic
      case ('dirichlet')
         type_id = bc_dirichlet
      case ('neumann', 'zero_gradient')
         type_id = bc_neumann
      case default
         call fatal_error('bc', 'unknown boundary type '//trim(text))
      end select
   end function parse_bc_type

   integer function patch_type_for_face(mesh, bc, face_id) result(type_id)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id

      integer :: patch_id

      patch_id = mesh%faces(face_id)%patch
      if (patch_id <= 0) then
         type_id = bc_unknown
      else
         type_id = bc%patches(patch_id)%type_id
      end if
   end function patch_type_for_face


   logical function is_periodic_face(mesh, bc, face_id)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id

      is_periodic_face = patch_type_for_face(mesh, bc, face_id) == bc_periodic
   end function is_periodic_face


   integer function face_effective_neighbor(mesh, bc, face_id, cell_id) result(neighbor)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      integer, intent(in) :: cell_id

      if (mesh%faces(face_id)%owner == cell_id) then
         neighbor = mesh%faces(face_id)%neighbor
      else
         neighbor = mesh%faces(face_id)%owner
      end if

      if (neighbor == 0 .and. is_periodic_face(mesh, bc, face_id)) then
         neighbor = mesh%faces(face_id)%periodic_neighbor
      end if
   end function face_effective_neighbor


   !> Evaluate velocity at a boundary face.
   !!
   !! @param mesh Mesh data structure.
   !! @param bc Boundary condition set.
   !! @param face_id ID of the boundary face.
   !! @param interior_velocity Velocity in the owner cell.
   !! @param value Resulting velocity at the face.
   subroutine boundary_velocity(mesh, bc, face_id, interior_velocity, value)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      real(rk), intent(in) :: interior_velocity(3)
      real(rk), intent(out) :: value(3)

      integer :: patch_id
      real(rk) :: nhat(3)
      real(rk) :: un

      patch_id = mesh%faces(face_id)%patch
      if (patch_id <= 0) then
         value = interior_velocity
         return
      end if

      select case (bc%patches(patch_id)%type_id)
      case (bc_wall, bc_dirichlet)
         value = bc%patches(patch_id)%velocity
      case (bc_symmetry)
         nhat = mesh%faces(face_id)%normal
         un = dot_product(interior_velocity, nhat)
         value = interior_velocity - un * nhat
      case (bc_neumann)
         value = interior_velocity
      case (bc_periodic)
         value = interior_velocity
      case default
         value = interior_velocity
      end select
   end subroutine boundary_velocity


   subroutine ensure_periodic_links(mesh, patch_id)
      type(mesh_t), intent(in) :: mesh
      integer, intent(in) :: patch_id

      integer :: i, face_id

      do i = 1, mesh%patches(patch_id)%nfaces
         face_id = mesh%patches(patch_id)%face_ids(i)
         if (mesh%faces(face_id)%periodic_neighbor <= 0) then
            call fatal_error('bc', 'periodic patch '//trim(mesh%patches(patch_id)%name)// &
                             ' has no periodic.dat face links')
         end if
      end do
   end subroutine ensure_periodic_links

   !> Evaluate species mass fractions at a boundary face.
   !!
   !! @param mesh Mesh data structure.
   !! @param bc Boundary condition set.
   !! @param face_id ID of the boundary face.
   !! @param k Species index.
   !! @param interior_Y Mass fraction in the owner cell.
   !! @param ext_Y Resulting mass fraction at the face (or external ghost state).
   !! @param is_dirichlet Flag indicating if the boundary is Dirichlet.
   subroutine boundary_species(mesh, bc, face_id, k, interior_Y, ext_Y, is_dirichlet)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      integer, intent(in) :: k
      real(rk), intent(in) :: interior_Y
      real(rk), intent(out) :: ext_Y
      logical, intent(out) :: is_dirichlet

      integer :: patch_id

      patch_id = mesh%faces(face_id)%patch
      if (patch_id <= 0) then
         ext_Y = interior_Y
         is_dirichlet = .false.
         return
      end if

      select case (bc%patches(patch_id)%species_type_id)
      case (bc_dirichlet)
         ext_Y = bc%patches(patch_id)%species_Y(k)
         is_dirichlet = .true.
      case default
         ext_Y = interior_Y
         is_dirichlet = .false.
      end select
   end subroutine boundary_species

end module mod_bc

