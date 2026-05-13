!> Finite-Volume transport solver for chemical species mass fractions.
!!
!! This module implements the solution of the transport equation for \(N\)
!! chemical species mass fractions \(Y_k\). The solver supports:
!! 1. **Upwind Advection**: A 1st-order stable scheme for robust transport 
!!    of sharp scalar gradients.
!! 2. **Corrected Diffusion**: Diffusive fluxes are explicitly corrected 
!!    to ensure the net mass flux sums to zero \(\sum \mathbf{j}_k = 0\).
!! 3. **Mass Conservation**: Enforces \(\sum Y_k = 1\) and boundedness \([0, 1]\)
!!    after every timestep.
!! 4. **MPI Synchronization**: Efficiently gathers owned-cell updates into 
!!    the globally replicated mesh field.
module mod_species
   use mod_kinds, only : rk, zero, one, fatal_error, name_len, lowercase
   use mod_mesh_types, only : mesh_t
   use mod_input, only : case_params_t
   use mod_bc, only : bc_set_t, bc_periodic, patch_type_for_face, face_effective_neighbor, boundary_species
   use mod_fields, only : flow_fields_t
   use mod_mpi_flow, only : flow_mpi_t, flow_exchange_cell_matrix
   use mod_flow_projection, only : face_normal_distance
   use mod_transport_properties, only : transport_properties_t
   implicit none

   private

   public :: species_fields_t
   public :: initialize_species, finalize_species, advance_species_transport

   !> Container for multi-species mass fraction fields.
   type :: species_fields_t
      integer :: nspecies = 0                       !< Total number of transport species \(N_s\).
      real(rk), allocatable :: Y(:,:)               !< Current mass fractions \(Y_k\) (nspecies, ncells).
      real(rk), allocatable :: Y_old(:,:)           !< Mass fractions from previous step \(n\).
      character(len=name_len), allocatable :: names(:) !< Array of species names (e.g., "H2", "O2", "N2").
   end type species_fields_t

contains

   !> Populates species fields with initial mass fractions and handles naming.
   !!
   !! Performs name-based matching between the transport registry and 
   !! namelist initial conditions. Normalizes the initial mixture to 
   !! ensure the physical constraint \(\sum Y_k = 1\) is met at \(t=0\).
   !!
   !! @param mesh The computational mesh.
   !! @param params Input configuration.
   !! @param species The fields to initialize.
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

      ! Match namelist initial conditions by species name
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

      ! Normalize the initial mixture vector
      sum_Y = sum(init_mixture)
      if (sum_Y > zero) then
         init_mixture = init_mixture / sum_Y
      else
         ! Fallback: If no IC specified, set the mixture to 100% of the first species.
         if (species%nspecies > 0) init_mixture(1) = one
      end if

      do c = 1, mesh%ncells
         species%Y(:, c) = init_mixture
      end do

      species%Y_old = species%Y
   end subroutine initialize_species


   !> Safely deallocates species fields and names.
   subroutine finalize_species(species)
      type(species_fields_t), intent(inout) :: species

      if (allocated(species%Y)) deallocate(species%Y)
      if (allocated(species%Y_old)) deallocate(species%Y_old)
      if (allocated(species%names)) deallocate(species%names)
      species%nspecies = 0
   end subroutine finalize_species


   !> Performs one explicit Euler step for species transport.
   !!
   !! This routine evaluates:
   !! - **Advective Fluxes**: \(J_{adv,k} = \dot{m}_f Y_{upwind,k}\).
   !! - **Diffusive Fluxes**: \(J_{diff,k} = -D_k A_f \nabla Y_k\).
   !! - **Flux Correction**: Subtracts \(Y_k \sum J_{diff,k}\) to satisfy mass conservation.
   !! - **Bounding & Normalization**: Clamps \(Y_k \in [0, 1]\) and enforces \(\sum Y_k = 1\) locally.
   !!
   !! @param mesh The computational mesh.
   !! @param flow MPI decomposition data for synchronization.
   !! @param bc Boundary condition settings.
   !! @param params Simulation parameters (dt, etc).
   !! @param fields Flow field (velocity/face fluxes).
   !! @param species Mass fraction fields to update.
   !! @param transport Physical properties (diffusivities).
   subroutine advance_species_transport(mesh, flow, bc, params, fields, species, transport)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(inout) :: flow
      type(bc_set_t), intent(in) :: bc
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(species_fields_t), intent(inout) :: species
      type(transport_properties_t), intent(in) :: transport

      real(rk), allocatable :: dY(:), diff_flux(:), adv_flux(:), Y_face_lin(:)
      real(rk) :: flux, face_area, dist, sum_diff_flux
      real(rk) :: Y_cell, Y_other, Y_face
      real(rk) :: diff_face
      integer :: c, f, fid, neighbor
      integer :: k
      real(rk) :: sum_Y
      logical :: is_dirichlet

      if (species%nspecies <= 0) return

      allocate(dY(species%nspecies))
      allocate(diff_flux(species%nspecies))
      allocate(adv_flux(species%nspecies))
      allocate(Y_face_lin(species%nspecies))

      species%Y_old = species%Y

      ! Iterate through MPI-owned cells
      do c = flow%first_cell, flow%last_cell
         dY = zero

         do f = 1, mesh%ncell_faces(c)
            fid = mesh%cell_faces(f,c)
            if (mesh%faces(fid)%owner == c) then
               flux = fields%face_flux(fid)
            else
               flux = -fields%face_flux(fid)
            end if
            neighbor = face_effective_neighbor(mesh, bc, fid, c)

            face_area = mesh%faces(fid)%area
            dist = face_normal_distance(mesh, bc, fid, c, neighbor)

            sum_diff_flux = zero
            do k = 1, species%nspecies
               Y_cell = species%Y_old(k, c)
               
               if (neighbor == 0) then
                  call boundary_species(mesh, bc, fid, k, Y_cell, Y_other, is_dirichlet)
               else
                  Y_other = species%Y_old(k, neighbor)
                  is_dirichlet = .true.
               end if

               ! 1. Advective flux using Upwind discretization.
               ! flux is oriented outward from the current cell.
               if (flux > zero) then
                  Y_face = Y_cell
               else
                  Y_face = Y_other
               end if
               adv_flux(k) = -flux * Y_face

               ! 2. Diffusive flux using central difference
               diff_flux(k) = zero
               if (neighbor /= 0 .or. is_dirichlet) then
                  if (neighbor == 0) then
                     diff_face = transport%diffusivity(k, c)
                  else
                     diff_face = 0.5_rk * (transport%diffusivity(k, c) + transport%diffusivity(k, neighbor))
                  end if
                  diff_flux(k) = diff_face * (Y_other - Y_cell) / dist * face_area
               end if
               
               sum_diff_flux = sum_diff_flux + diff_flux(k)
               Y_face_lin(k) = 0.5_rk * (Y_cell + Y_other)
            end do

            ! 3. Apply Correction Velocity to ensure mass conservation: sum(j_k) = 0
            do k = 1, species%nspecies
               dY(k) = dY(k) + adv_flux(k) + (diff_flux(k) - Y_face_lin(k) * sum_diff_flux)
            end do
         end do

         ! Explicit timestep update
         do k = 1, species%nspecies
            species%Y(k,c) = species%Y_old(k,c) + params%dt * dY(k) / mesh%cells(c)%volume
            
            ! Ensure boundedness: Y_k must remain in [0, 1]
            if (species%Y(k,c) < zero) species%Y(k,c) = zero
            if (species%Y(k,c) > one)  species%Y(k,c) = one
         end do

         ! Local renormalization: sum(Y_k) = 1
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
      deallocate(diff_flux)
      deallocate(adv_flux)
      deallocate(Y_face_lin)

      ! Synchronize updated species ghosts for the next transport/property step.
      call flow_exchange_cell_matrix(flow, species%Y)

   end subroutine advance_species_transport


end module mod_species
