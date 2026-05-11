!> Input parameter parsing and validation.
!!
!! This module handles reading the case.nml namelist file and storing
!! all physical and solver-related parameters.
module mod_input
   use mod_kinds, only : rk, zero, one, name_len, path_len, fatal_error, lowercase
   implicit none

   private

   integer, parameter, public :: max_patches = 64
   integer, parameter, public :: max_species = 16

   !> Container for all case parameters.
   type, public :: case_params_t
      character(len=path_len) :: mesh_dir = "mesh_native" !< Directory containing mesh .dat files

      integer :: nsteps = 100            !< Number of timesteps to run
      real(rk) :: dt = 1.0e-3_rk         !< Timestep size [s]
      integer :: output_interval = 10    !< Steps between VTU/diagnostic outputs

      real(rk) :: rho = one              !< Constant density [kg/m^3]
      real(rk) :: nu = 1.0e-2_rk         !< Constant kinematic viscosity [m^2/s]

      integer :: pressure_max_iter = 300 !< Maximum CG iterations for pressure Poisson solve
      real(rk) :: pressure_tol = 1.0e-10_rk !< Tolerance for pressure Poisson solve
      real(rk) :: body_force(3) = zero   !< Volumetric body force [N/m^3]
      character(len=name_len) :: convection_scheme = "upwind" !< Advection scheme ("upwind" or "central")

      integer :: n_patches = 0           !< Number of boundary patches
      character(len=name_len) :: patch_name(max_patches) = "" !< Names of patches
      character(len=name_len) :: patch_type(max_patches) = "" !< Legacy patch types

      !> Stage 2 boundary condition split.
      !! These are optional in case.nml. If left blank, mod_bc can fall back
      !! to the legacy patch_type behavior.
      character(len=name_len) :: patch_velocity_type(max_patches) = "" !< Velocity BC type
      character(len=name_len) :: patch_pressure_type(max_patches) = "" !< Pressure BC type

      real(rk) :: patch_u(max_patches) = zero    !< Specified u-velocity on patch
      real(rk) :: patch_v(max_patches) = zero    !< Specified v-velocity on patch
      real(rk) :: patch_w(max_patches) = zero    !< Specified w-velocity on patch
      real(rk) :: patch_p(max_patches) = zero    !< Specified pressure on patch
      real(rk) :: patch_dpdn(max_patches) = zero !< Specified pressure gradient on patch

      !> Stage 3B species boundary conditions
      character(len=name_len) :: patch_species_type(max_patches) = "" !< Species BC type
      real(rk) :: patch_Y(max_species, max_patches) = zero           !< Specified mass fractions on patch


      character(len=path_len) :: output_dir = "output" !< Directory for output files
      logical :: write_vtu = .true.                    !< Enable VTU output
      logical :: write_diagnostics = .true.            !< Enable diagnostics.csv output

      ! Species parameters
      integer :: nspecies = 0                                     !< Number of transport species
      character(len=name_len) :: species_name(max_species) = ""   !< Names of species
      real(rk) :: species_diffusivity(max_species) = 0.0_rk       !< Constant diffusivity if Cantera disabled
      real(rk) :: initial_Y(max_species) = 0.0_rk                 !< Initial mass fractions in domain
      logical :: enable_cantera = .false.                         !< Enable Cantera property bridge
      character(len=path_len) :: cantera_mech_file = "gri30.yaml" !< Cantera mechanism file
      real(rk) :: background_temp = 300.0_rk                      !< Background temperature for properties [K]
      real(rk) :: background_press = 101325.0_rk                  !< Background pressure for properties [Pa]
   end type case_params_t

   public :: read_case_params

contains

   !> Read all parameters from a namelist file.
   !!
   !! @param filename Path to the case.nml file.
   !! @param params Parameters structure to populate.
   subroutine read_case_params(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      call read_mesh_input(filename, params)
      call read_time_input(filename, params)
      call read_fluid_input(filename, params)
      call read_solver_input(filename, params)
      call read_boundary_input(filename, params)
      call read_species_input(filename, params)
      call read_output_input(filename, params)
      call validate_params(params)
   end subroutine read_case_params


   subroutine read_mesh_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      character(len=path_len) :: mesh_dir
      integer :: unit_id, ios

      namelist /mesh_input/ mesh_dir

      mesh_dir = params%mesh_dir

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=mesh_input, iostat=ios)
         close(unit_id)
      end if

      if (ios > 0) call fatal_error('input', 'failed reading &mesh_input')

      params%mesh_dir = mesh_dir
   end subroutine read_mesh_input


   subroutine read_time_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      integer :: nsteps, output_interval
      real(rk) :: dt
      integer :: unit_id, ios

      namelist /time_input/ nsteps, dt, output_interval

      nsteps = params%nsteps
      dt = params%dt
      output_interval = params%output_interval

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=time_input, iostat=ios)
         close(unit_id)
      end if

      if (ios > 0) call fatal_error('input', 'failed reading &time_input')

      params%nsteps = nsteps
      params%dt = dt
      params%output_interval = output_interval
   end subroutine read_time_input


   subroutine read_fluid_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      real(rk) :: rho, nu
      integer :: unit_id, ios

      namelist /fluid_input/ rho, nu

      rho = params%rho
      nu = params%nu

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=fluid_input, iostat=ios)
         close(unit_id)
      end if

      if (ios > 0) call fatal_error('input', 'failed reading &fluid_input')

      params%rho = rho
      params%nu = nu
   end subroutine read_fluid_input


   subroutine read_solver_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      integer :: pressure_max_iter
      real(rk) :: pressure_tol
      real(rk) :: body_force_x, body_force_y, body_force_z
      character(len=name_len) :: convection_scheme
      integer :: unit_id, ios

      namelist /solver_input/ pressure_max_iter, pressure_tol, &
                              body_force_x, body_force_y, body_force_z, &
                              convection_scheme

      pressure_max_iter = params%pressure_max_iter
      pressure_tol = params%pressure_tol
      body_force_x = params%body_force(1)
      body_force_y = params%body_force(2)
      body_force_z = params%body_force(3)
      convection_scheme = params%convection_scheme

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=solver_input, iostat=ios)
         close(unit_id)
      end if

      if (ios > 0) call fatal_error('input', 'failed reading &solver_input')

      params%pressure_max_iter = pressure_max_iter
      params%pressure_tol = pressure_tol
      params%body_force = [body_force_x, body_force_y, body_force_z]
      params%convection_scheme = trim(lowercase(convection_scheme))
   end subroutine read_solver_input


   subroutine read_boundary_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      integer :: n_patches
      character(len=name_len) :: patch_name(max_patches)
      character(len=name_len) :: patch_type(max_patches)
      character(len=name_len) :: patch_velocity_type(max_patches)
      character(len=name_len) :: patch_pressure_type(max_patches)
      character(len=name_len) :: patch_species_type(max_patches)
      real(rk) :: patch_u(max_patches)
      real(rk) :: patch_v(max_patches)
      real(rk) :: patch_w(max_patches)
      real(rk) :: patch_p(max_patches)
      real(rk) :: patch_dpdn(max_patches)
      real(rk) :: patch_Y(max_species, max_patches)
      integer :: unit_id, ios, i

      namelist /boundary_input/ n_patches, patch_name, patch_type, &
                                 patch_velocity_type, patch_pressure_type, &
                                 patch_species_type, &
                                 patch_u, patch_v, patch_w, patch_p, patch_dpdn, patch_Y

      n_patches = params%n_patches
      patch_name = params%patch_name
      patch_type = params%patch_type
      patch_velocity_type = params%patch_velocity_type
      patch_pressure_type = params%patch_pressure_type
      patch_species_type = params%patch_species_type
      patch_u = params%patch_u
      patch_v = params%patch_v
      patch_w = params%patch_w
      patch_p = params%patch_p
      patch_dpdn = params%patch_dpdn
      patch_Y = params%patch_Y

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=boundary_input, iostat=ios)
         close(unit_id)
      end if

      if (ios > 0) call fatal_error('input', 'failed reading &boundary_input')

      params%n_patches = n_patches
      params%patch_name = patch_name
      params%patch_type = patch_type
      params%patch_velocity_type = patch_velocity_type
      params%patch_pressure_type = patch_pressure_type
      params%patch_species_type = patch_species_type

      do i = 1, max_patches
         params%patch_type(i) = trim(lowercase(params%patch_type(i)))
         params%patch_velocity_type(i) = trim(lowercase(params%patch_velocity_type(i)))
         params%patch_pressure_type(i) = trim(lowercase(params%patch_pressure_type(i)))
         params%patch_species_type(i) = trim(lowercase(params%patch_species_type(i)))
      end do

      params%patch_u = patch_u
      params%patch_v = patch_v
      params%patch_w = patch_w
      params%patch_p = patch_p
      params%patch_dpdn = patch_dpdn
      params%patch_Y = patch_Y
   end subroutine read_boundary_input


   subroutine read_output_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      character(len=path_len) :: output_dir
      logical :: write_vtu, write_diagnostics
      integer :: unit_id, ios

      namelist /output_input/ output_dir, write_vtu, write_diagnostics

      output_dir = params%output_dir
      write_vtu = params%write_vtu
      write_diagnostics = params%write_diagnostics

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=output_input, iostat=ios)
         close(unit_id)
      end if

      if (ios > 0) call fatal_error('input', 'failed reading &output_input')

      params%output_dir = output_dir
      params%write_vtu = write_vtu
      params%write_diagnostics = write_diagnostics
   end subroutine read_output_input


   subroutine read_species_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      integer :: nspecies
      character(len=name_len) :: species_name(max_species)
      real(rk) :: species_diffusivity(max_species)
      real(rk) :: initial_Y(max_species)
      logical :: enable_cantera
      character(len=path_len) :: cantera_mech_file
      real(rk) :: background_temp
      real(rk) :: background_press
      integer :: unit_id, ios

      namelist /species_input/ nspecies, species_name, species_diffusivity, initial_Y, &
                               enable_cantera, cantera_mech_file, background_temp, background_press

      nspecies = params%nspecies
      species_name = params%species_name
      species_diffusivity = params%species_diffusivity
      initial_Y = params%initial_Y
      enable_cantera = params%enable_cantera
      cantera_mech_file = params%cantera_mech_file
      background_temp = params%background_temp
      background_press = params%background_press

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=species_input, iostat=ios)
         close(unit_id)
      end if

      ! If ios > 0, it might just mean the &species_input block doesn't exist.
      ! We can treat species as optional.
      if (ios == 0) then
         params%nspecies = nspecies
         params%species_name = species_name
         params%species_diffusivity = species_diffusivity
         params%initial_Y = initial_Y
         params%enable_cantera = enable_cantera
         params%cantera_mech_file = cantera_mech_file
         params%background_temp = background_temp
         params%background_press = background_press
      end if
      if (ios > 0) then
         print *, "read_species_input failed, ios = ", ios
      end if
   end subroutine read_species_input


   subroutine open_namelist_file(filename, unit_id, ios)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: unit_id
      integer, intent(out) :: ios

      open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)

      if (ios /= 0) unit_id = -1
   end subroutine open_namelist_file


   subroutine validate_params(params)
      type(case_params_t), intent(in) :: params

      if (params%nsteps < 0) call fatal_error('input', 'nsteps must be non-negative')
      if (params%dt <= zero) call fatal_error('input', 'dt must be positive')
      if (params%output_interval <= 0) call fatal_error('input', 'output_interval must be positive')
      if (params%rho <= zero) call fatal_error('input', 'rho must be positive')
      if (params%nu < zero) call fatal_error('input', 'nu must be non-negative')
      if (params%pressure_max_iter <= 0) call fatal_error('input', 'pressure_max_iter must be positive')
      if (params%pressure_tol <= zero) call fatal_error('input', 'pressure_tol must be positive')

      if (params%n_patches < 0 .or. params%n_patches > max_patches) then
         call fatal_error('input', 'n_patches is outside supported range')
      end if

      if (len_trim(params%mesh_dir) == 0) then
         call fatal_error('input', 'mesh_dir cannot be empty')
      end if

      if (len_trim(params%output_dir) == 0) then
         call fatal_error('input', 'output_dir cannot be empty')
      end if

      if (params%nspecies < 0 .or. params%nspecies > max_species) then
         call fatal_error('input', 'nspecies is outside supported range')
      end if

      call validate_boundary_arrays(params)
   end subroutine validate_params


   subroutine validate_boundary_arrays(params)
      type(case_params_t), intent(in) :: params

      integer :: i

      do i = 1, params%n_patches
         if (len_trim(params%patch_name(i)) == 0) then
            call fatal_error('input', 'patch_name entry cannot be empty')
         end if

         if (len_trim(params%patch_type(i)) == 0) then
            call fatal_error('input', 'patch_type entry cannot be empty')
         end if
      end do
   end subroutine validate_boundary_arrays

end module mod_input