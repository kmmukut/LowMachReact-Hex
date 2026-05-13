!> Boundary condition (BC) management, parsing, and physical evaluation.
!!
!! This module provides the infrastructure to map user-defined case 
!! configuration (`case.nml`) to the geometric patches of the mesh. 
!! It handles the storage of boundary data (velocity, pressure, species) 
!! and provides low-level evaluation routines used by the spatial operators 
!! (e.g., flux calculation, gradient evaluation).
!!
!! Supported BC Types:
!! 1. **`bc_wall`**: No-slip/No-penetration condition. 
!!    - Velocity is set to the specified wall velocity (default zero).
!!    - Pressure gradient is zero (Neumann).
!! 2. **`bc_symmetry`**: Zero-gradient for parallel components, zero for normal.
!!    - Enforces no-flow across the boundary while allowing slip.
!! 3. **`bc_periodic`**: Topology-linked faces for repeating domains.
!!    - Requires a valid `periodic.dat` file created by the mesh converter.
!! 4. **`bc_dirichlet`**: Fixed value (Inlet).
!!    - Explicitly sets the field value at the face.
!! 5. **`bc_neumann`**: Fixed gradient (Zero-Gradient Outlet).
!!    - Sets the face value equal to the adjacent interior cell value.
module mod_bc
   use mod_kinds, only : rk, zero, name_len, fatal_error, lowercase
   use mod_input, only : case_params_t, max_species
   use mod_mesh_types, only : mesh_t
   implicit none

   private

   integer, parameter, public :: bc_unknown   = 0 !< Unspecified or invalid BC type.
   integer, parameter, public :: bc_wall      = 1 !< Standard no-slip solid wall.
   integer, parameter, public :: bc_symmetry  = 2 !< Symmetry plane / slip wall.
   integer, parameter, public :: bc_periodic  = 3 !< Periodic (cyclic) boundary.
   integer, parameter, public :: bc_dirichlet = 4 !< Fixed value (e.g., Inlet).
   integer, parameter, public :: bc_neumann   = 5 !< Fixed gradient (e.g., Zero-Gradient Outlet).

   !> Container for boundary data assigned to a specific mesh patch.
   !!
   !! Each patch in the mesh can have different BC types for different 
   !! physical fields (U, P, Y). This structure stores the numeric type 
   !! IDs and the associated physical values.
   type, public :: bc_patch_t
      integer :: patch_id = 0                    !< Link to the geometric `mesh%patch(id)`.
      character(len=name_len) :: name = ""       !< Human-readable name (e.g., "inlet").
      character(len=name_len) :: type_name = ""  !< Input type string from namelist (e.g., "wall").
      integer :: type_id = bc_unknown            !< Master BC type ID for the patch.
      real(rk) :: velocity(3) = zero             !< Specified velocity vector \((u,v,w)\) [m/s].
      real(rk) :: pressure = zero                !< Specified static pressure \(P\) [Pa].
      real(rk) :: dpdn = zero                    !< Specified pressure gradient \(dP/dn\) [Pa/m].

      !> Species boundary settings.
      integer :: species_type_id = bc_unknown    !< BC type applied to species transport.
      real(rk) :: species_Y(max_species) = zero  !< Specified mass fractions \(Y_k\) for Dirichlet boundaries.

      !> Field-specific overrides.
      integer :: velocity_type_id = bc_unknown   !< BC type for velocity (if different from master).
      integer :: pressure_type_id = bc_unknown   !< BC type for pressure (if different from master).

      !> Temperature/enthalpy boundary settings.
      integer :: temperature_type_id = bc_unknown !< BC type applied to temperature/enthalpy transport.
      real(rk) :: temperature = 300.0_rk          !< Specified boundary temperature [K].
   end type bc_patch_t

   !> Global set of boundary conditions covering all mesh patches.
   type :: bc_set_t
      integer :: npatches = 0                        !< Total number of patches defined in the mesh.
      type(bc_patch_t), allocatable :: patches(:)    !< Array of boundary settings for each patch.
   end type bc_set_t

   public :: build_bc_set, finalize_bc_set, bc_set_t
   public :: patch_type_for_face, boundary_velocity, face_effective_neighbor
   public :: boundary_pressure_type, boundary_pressure
   public :: boundary_species
   public :: boundary_temperature
   public :: is_periodic_face

contains

   !> Synchronizes namelist parameters with mesh patches to create a complete BC set.
   !!
   !! This routine iterates through all patches in the `mesh_t` structure and 
   !! searches for matching names in the `case_params_t` object. It converts 
   !! string-based inputs into internal type IDs.
   !!
   !! @param mesh The computational mesh containing patch definitions.
   !! @param params Parsed case configuration from `case.nml`.
   !! @param bc The boundary condition set to be populated.
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
               bc%patches(p)%temperature = params%patch_T(q)

               if (trim(params%patch_velocity_type(q)) /= "") then
                  bc%patches(p)%velocity_type_id = parse_bc_type(params%patch_velocity_type(q))
               else
                  bc%patches(p)%velocity_type_id = parse_bc_type(params%patch_type(q))
               end if

               if (trim(params%patch_pressure_type(q)) /= "") then
                  bc%patches(p)%pressure_type_id = parse_bc_type(params%patch_pressure_type(q))
               else
                  bc%patches(p)%pressure_type_id = parse_bc_type(params%patch_type(q))
               end if

               if (trim(params%patch_species_type(q)) /= "") then
                  bc%patches(p)%species_type_id = parse_bc_type(params%patch_species_type(q))
               else
                  bc%patches(p)%species_type_id = parse_bc_type(params%patch_type(q))
               end if
               bc%patches(p)%species_Y(:) = params%patch_Y(:, q)

               if (trim(params%patch_temperature_type(q)) /= "") then
                  bc%patches(p)%temperature_type_id = parse_bc_type(params%patch_temperature_type(q))
               else
                  bc%patches(p)%temperature_type_id = parse_bc_type(params%patch_type(q))
               end if

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


   !> Safely deallocates the boundary condition patch array.
   subroutine finalize_bc_set(bc)
      type(bc_set_t), intent(inout) :: bc
      if (allocated(bc%patches)) deallocate(bc%patches)
      bc%npatches = 0
   end subroutine finalize_bc_set


   !> Converts a case-insensitive string to its corresponding internal BC type ID.
   !!
   !! @param text String representation (e.g., "Wall", "Periodic", "Dirichlet").
   !!       Supports legacy aliases: "fixed_value" (Dirichlet), "zero_gradient" (Neumann), 
   !!       "no_slip" (Wall), "slip" (Symmetry).
   !! @result The integer type ID (e.g., `bc_wall`).
   integer function parse_bc_type(text) result(type_id)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: key

      type_id = bc_unknown

      key = trim(lowercase(text))
      select case (trim(key))
      case ('wall', 'no_slip', 'moving_wall')
         type_id = bc_wall
      case ('symmetry', 'symmetric', 'slip')
         type_id = bc_symmetry
      case ('periodic')
         type_id = bc_periodic
      case ('dirichlet', 'fixed_value')
         type_id = bc_dirichlet
      case ('neumann', 'zero_gradient')
         type_id = bc_neumann
      case default
         call fatal_error('bc', 'unknown boundary type '//trim(text))
      end select
   end function parse_bc_type

   !> Retrieves the master BC type for a given face ID.
   !!
   !! @param mesh The computational mesh.
   !! @param bc The active BC set.
   !! @param face_id Global index of the face.
   !! @result The BC type ID of the patch containing the face.
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


   !> Returns true if the face belongs to a periodic boundary.
   logical function is_periodic_face(mesh, bc, face_id)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      is_periodic_face = patch_type_for_face(mesh, bc, face_id) == bc_periodic
   end function is_periodic_face


   !> Returns the neighbor cell index, accounting for periodic connectivity.
   !!
   !! If the face is on a periodic boundary, this returns the `periodic_neighbor` 
   !! stored in the face structure. Otherwise, it returns the standard neighbor.
   !!
   !! @param mesh The computational mesh.
   !! @param bc The active BC set.
   !! @param face_id Global index of the face.
   !! @param cell_id Index of the cell requesting the neighbor.
   !! @result The effective neighbor cell index.
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


   !> Evaluates the velocity vector at a boundary face.
   !!
   !! Implements standard FV boundary evaluations:
   !! - **Wall/Dirichlet**: Returns the specified patch velocity.
   !! - **Symmetry**: Subtracts the normal component from the interior velocity.
   !! - **Neumann/Periodic**: Extrapolates the interior velocity to the face.
   !!
   !! @param mesh Mesh data structure.
   !! @param bc Boundary condition set.
   !! @param face_id ID of the boundary face.
   !! @param interior_velocity Velocity in the owner cell [m/s].
   !! @param value Resulting velocity at the face [m/s].
   subroutine boundary_velocity(mesh, bc, face_id, interior_velocity, value)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      real(rk), intent(in) :: interior_velocity(3)
      real(rk), intent(out) :: value(3)

      integer :: patch_id
      real(rk) :: nhat(3), un

      patch_id = mesh%faces(face_id)%patch
      if (patch_id <= 0) then
         value = interior_velocity
         return
      end if

      select case (bc%patches(patch_id)%velocity_type_id)
      case (bc_wall, bc_dirichlet)
         value = bc%patches(patch_id)%velocity
      case (bc_symmetry)
         nhat = mesh%faces(face_id)%normal
         un = dot_product(interior_velocity, nhat)
         value = interior_velocity - un * nhat
      case default
         value = interior_velocity
      end select
   end subroutine boundary_velocity


   !> Validates that periodic patches have correctly established links.
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

   !> Evaluates species mass fractions at a boundary face.
   !!
   !! Currently supports fixed-value (Dirichlet) and zero-gradient (Neumann) BCs.
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

   !> Returns the pressure BC type for a given face.

   !> Evaluates temperature at a boundary face.
   !!
   !! Dirichlet/fixed_value boundaries return the configured patch_T.
   !! All other boundary types currently behave as zero-gradient boundaries.
   subroutine boundary_temperature(mesh, bc, face_id, interior_T, ext_T, is_dirichlet)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      real(rk), intent(in) :: interior_T
      real(rk), intent(out) :: ext_T
      logical, intent(out) :: is_dirichlet

      integer :: patch_id

      patch_id = mesh%faces(face_id)%patch
      if (patch_id <= 0) then
         ext_T = interior_T
         is_dirichlet = .false.
         return
      end if

      select case (bc%patches(patch_id)%temperature_type_id)
      case (bc_dirichlet)
         ext_T = bc%patches(patch_id)%temperature
         is_dirichlet = .true.
      case default
         ext_T = interior_T
         is_dirichlet = .false.
      end select
   end subroutine boundary_temperature

   integer function boundary_pressure_type(mesh, bc, face_id) result(type_id)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id

      integer :: patch_id

      patch_id = mesh%faces(face_id)%patch
      if (patch_id <= 0) then
         type_id = bc_unknown
      else
         type_id = bc%patches(patch_id)%pressure_type_id
      end if
   end function boundary_pressure_type

   !> Evaluates pressure at a boundary face.
   subroutine boundary_pressure(mesh, bc, face_id, interior_p, ext_p, is_dirichlet)
      type(mesh_t), intent(in) :: mesh
      type(bc_set_t), intent(in) :: bc
      integer, intent(in) :: face_id
      real(rk), intent(in) :: interior_p
      real(rk), intent(out) :: ext_p
      logical, intent(out) :: is_dirichlet

      integer :: patch_id

      patch_id = mesh%faces(face_id)%patch
      if (patch_id <= 0) then
         ext_p = interior_p
         is_dirichlet = .false.
         return
      end if

      select case (bc%patches(patch_id)%pressure_type_id)
      case (bc_dirichlet)
         ext_p = bc%patches(patch_id)%pressure
         is_dirichlet = .true.
      case default
         ext_p = interior_p
         is_dirichlet = .false.
      end select
   end subroutine boundary_pressure

end module mod_bc

