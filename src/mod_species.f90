!> Species transport solver for passive scalars.
!!
!! This module implements the finite-volume transport of multiple species
!! mass fractions using explicit time integration and upwind advection.
module mod_species
   use mod_kinds, only : rk, zero, one, fatal_error, name_len, lowercase
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t
   use mod_bc, only : bc_set_t, bc_periodic, patch_type_for_face, face_effective_neighbor, boundary_species
   use mod_fields, only : flow_fields_t
   use mod_mpi_flow, only : flow_mpi_t, flow_allgather_owned_scalar
   use mod_flow_projection, only : face_normal_distance
   use mod_transport_properties, only : transport_properties_t
   implicit none

   private

   public :: species_fields_t
   public :: initialize_species, finalize_species, advance_species_transport

   !> Data structure for cell-centered species mass fractions.
   type :: species_fields_t
      integer :: nspecies = 0                       !< Number of species
      real(rk), allocatable :: Y(:,:)               !< Current mass fractions (nspecies, ncells)
      real(rk), allocatable :: Y_old(:,:)           !< Previous-step mass fractions (nspecies, ncells)
      character(len=name_len), allocatable :: names(:) !< Species names
   end type species_fields_t

contains

   !> Initialize species mass fractions from case parameters.
   !!
   !! @param mesh Full replicated mesh.
   !! @param params Case configuration parameters.
   !! @param species Species field data structure to initialize.
   subroutine initialize_species(mesh, params, species)
      type(mesh_t), intent(in) :: mesh
      type(case_params_t), intent(in) :: params
      type(species_fields_t), intent(inout) :: species

      integer :: c, k, j
      real(rk) :: sum_Y
      real(rk) :: init_mixture(params%nspecies)
      character(len=name_len) :: target_name

      call finalize_species(species)

      species%nspecies = params%nspecies
      if (species%nspecies <= 0) return

      allocate(species%Y(species%nspecies, mesh%ncells))
      allocate(species%Y_old(species%nspecies, mesh%ncells))
      allocate(species%names(species%nspecies))

      species%names = params%species_name(1:species%nspecies)

      ! ------------------------------------------------------------
      ! Name-based matching for initialization
      ! ------------------------------------------------------------
      init_mixture = zero
      
      do j = 1, params%namelist_nspecies
         target_name = trim(lowercase(params%namelist_species_name(j)))
         if (len_trim(target_name) == 0) cycle
         
         do k = 1, species%nspecies
            if (trim(lowercase(species%names(k))) == target_name) then
               init_mixture(k) = params%initial_Y(j)
               exit
            end if
         end do
      end do

      ! Normalize
      sum_Y = sum(init_mixture)
      if (sum_Y > zero) then
         init_mixture = init_mixture / sum_Y
      else
         ! Fallback if nothing matched or sum was zero
         if (species%nspecies > 0) init_mixture(1) = one
      end if

      do c = 1, mesh%ncells
         species%Y(:, c) = init_mixture
      end do

      species%Y_old = species%Y
   end subroutine initialize_species


   !> Deallocate all species-related arrays.
   !!
   !! @param species Species field data structure to finalize.
   subroutine finalize_species(species)
      type(species_fields_t), intent(inout) :: species

      if (allocated(species%Y)) deallocate(species%Y)
      if (allocated(species%Y_old)) deallocate(species%Y_old)
      if (allocated(species%names)) deallocate(species%names)
      species%nspecies = 0
   end subroutine finalize_species


   !> Advance the species transport equation by one timestep.
   !!
   !! Computes advective and diffusive fluxes, performs an explicit update,
   !! clamps values to [0, 1], and renormalizes sum(Y)=1.
   !!
   !! @param mesh Full replicated mesh.
   !! @param flow MPI flow decomposition data.
   !! @param bc Boundary condition set.
   !! @param params Case configuration parameters.
   !! @param fields Flow fields (velocity/fluxes) for advection.
   !! @param species Species field data to update.
   !! @param transport Transport properties (diffusivity).
   subroutine advance_species_transport(mesh, flow, bc, params, fields, species, transport)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(species_fields_t), intent(inout) :: species
      type(transport_properties_t), intent(in) :: transport

      real(rk), allocatable :: dY(:)
      real(rk) :: flux, face_area, dist
      real(rk) :: Y_owner, Y_neighbor, Y_face
      real(rk) :: diff_face
      integer :: c, f, fid, owner, neighbor
      integer :: k
      real(rk) :: sum_Y
      logical :: is_dirichlet

      if (species%nspecies <= 0) return

      allocate(dY(species%nspecies))

      species%Y_old = species%Y

      ! Update each cell
      do c = flow%first_cell, flow%last_cell
         dY = zero

         do f = 1, mesh%ncell_faces(c)
            fid = mesh%cell_faces(f,c)
            flux = fields%face_flux(fid)

            owner = mesh%faces(fid)%owner
            neighbor = face_effective_neighbor(mesh, bc, fid, c)

            face_area = mesh%faces(fid)%area
            dist = face_normal_distance(mesh, bc, fid, c, neighbor)

            do k = 1, species%nspecies
               Y_owner = species%Y_old(k, owner)
               
               if (neighbor == 0) then
                  call boundary_species(mesh, bc, fid, k, Y_owner, Y_neighbor, is_dirichlet)
               else
                  Y_neighbor = species%Y_old(k, neighbor)
                  is_dirichlet = .false.
               end if

               ! Upwind advection
               if (flux > zero) then
                  Y_face = Y_owner
               else
                  Y_face = Y_neighbor
               end if

               ! Advective flux (outward positive)
               dY(k) = dY(k) - flux * Y_face

               ! Diffusive flux (central difference, outward positive)
               ! flux_diff = - D * grad(Y) dot n * area
               if (neighbor /= 0 .or. is_dirichlet) then
                  if (neighbor == 0) then
                     diff_face = transport%diffusivity(k, owner)
                  else
                     diff_face = 0.5_rk * (transport%diffusivity(k, owner) + transport%diffusivity(k, neighbor))
                  end if
                  dY(k) = dY(k) + diff_face * (Y_neighbor - Y_owner) / dist * face_area
               end if
            end do
         end do

         ! Explicit update
         do k = 1, species%nspecies
            species%Y(k,c) = species%Y_old(k,c) + params%dt * dY(k) / mesh%cells(c)%volume
            
            ! Clamp to [0, 1]
            if (species%Y(k,c) < zero) species%Y(k,c) = zero
            if (species%Y(k,c) > one)  species%Y(k,c) = one
         end do

         ! Renormalize mass fractions so sum(Y_k) = 1
         sum_Y = zero
         do k = 1, species%nspecies
            sum_Y = sum_Y + species%Y(k,c)
         end do

         if (sum_Y > zero) then
            do k = 1, species%nspecies
               species%Y(k,c) = species%Y(k,c) / sum_Y
            end do
         end if
      end do

      deallocate(dY)

      ! We need to sync the updated species across all ranks since the full mesh is replicated.
      call flow_allgather_owned_species(flow, species)

   end subroutine advance_species_transport

   !> Synchronize species mass fractions across all MPI ranks.
   !!
   !! Since the full mesh is replicated on each rank, we must gather
   !! owned-cell updates into the global arrays.
   !!
   !! @param flow MPI flow decomposition data.
   !! @param species Species field data structure.
   subroutine flow_allgather_owned_species(flow, species)
      type(flow_mpi_t), intent(inout) :: flow
      type(species_fields_t), intent(inout) :: species
      real(rk), allocatable :: local_Y(:), global_Y(:)
      integer :: k
      
      allocate(local_Y(size(species%Y, 2)), global_Y(size(species%Y, 2)))
      do k = 1, species%nspecies
         local_Y = species%Y(k,:)
         call flow_allgather_owned_scalar(flow, local_Y, global_Y)
         species%Y(k,:) = global_Y
      end do
      deallocate(local_Y, global_Y)
   end subroutine flow_allgather_owned_species

end module mod_species
