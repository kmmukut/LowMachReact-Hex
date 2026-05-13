!> Parsing and validation of simulation parameters from `case.nml`.
!!
!! This module handles the ingestion of user-defined parameters using Fortran 
!! namelists. It is structured to read parameters in logical groups:
!! - **`mesh_input`**: Path to the native hexahedral grid files.
!! - **`time_input`**: Step counts, timestep size, and CFL control.
!! - **`fluid_input`**: Constant/Cantera properties (density, viscosity).
!! - **`solver_input`**: Linear solver tolerances and numerical schemes.
!! - **`boundary_input`**: Mapping of BC types and values to mesh patches.
!! - **`species_input`**: Multi-species transport and chemical mechanism selection.
!! - **`output_input`**: Control over VTK and diagnostic file generation.
!! - **`profiling_input`**: Control over profiling enable/disable and nested reporting.
module mod_input
   use mod_kinds, only : rk, zero, one, name_len, path_len, fatal_error, lowercase
   implicit none

   private

   integer, parameter, public :: max_patches = 64  !< Maximum number of boundary patches supported.
   integer, parameter, public :: max_species = 256 !< Maximum number of chemical species supported.

   !> Global container for all parsed simulation parameters.
   !!
   !! This structure acts as the source of truth for the entire solver. 
   !! It is populated once at startup and accessed as a read-only object 
   !! by the various solver modules.
   type, public :: case_params_t
      !> @name Mesh Configuration
      character(len=path_len) :: mesh_dir = "mesh_native" !< Directory containing `points.dat`, `cells.dat`, etc.

      !> @name Time Stepping
      integer :: nsteps = 0              !< Total number of timesteps to execute.
      real(rk) :: dt = zero              !< Fixed timestep size [s].
      integer :: output_interval = 1     !< Frequency of VTK/PVD output (in steps).
      logical :: use_dynamic_dt = .false.!< If true, `dt` scales to maintain target `max_cfl`.
      real(rk) :: max_cfl = 0.5_rk       !< Target maximum CFL number for stability.

      !> @name Fluid Properties
      real(rk) :: rho = one              !< Constant flow/projection density [kg/m^3].
      logical :: enable_variable_density = .false. !< Reserved: variable-density/low-Mach flow is not enabled yet.
      real(rk) :: nu = 1.0e-2_rk         !< Constant kinematic viscosity [m^2/s].
      logical :: enable_variable_nu = .false.    !< Allow Cantera-updated flow viscosity/nu; default keeps validation Re fixed.
      integer :: transport_update_interval = 1 !< Cantera transport-property update interval for mu/D_k only [steps].

      !> @name Linear Solver & Numerics
      integer :: pressure_max_iter = 300 !< Maximum Conjugate Gradient iterations for the Poisson solver.
      real(rk) :: pressure_tol = 1.0e-10_rk !< Relative residual tolerance for pressure convergence.
      real(rk) :: body_force(3) = zero   !< Constant volumetric acceleration vector \((a_x, a_y, a_z)\) [m/s^2].
      character(len=name_len) :: convection_scheme = "upwind" !< Scheme for advection: "upwind" (stable) or "central" (accurate).

      !> @name Boundary Condition Mapping
      integer :: n_patches = 0           !< Total number of patches defined in the namelist.
      character(len=name_len) :: patch_name(max_patches) = "" !< Names identifying mesh patches.
      character(len=name_len) :: patch_type(max_patches) = "" !< Legacy BC type string.

      !> Field-specific BC Overrides
      character(len=name_len) :: patch_velocity_type(max_patches) = "" !< Override BC type for velocity.
      character(len=name_len) :: patch_pressure_type(max_patches) = "" !< Override BC type for pressure.
      character(len=name_len) :: patch_temperature_type(max_patches) = "" !< Override BC type for temperature/enthalpy.

      real(rk) :: patch_u(max_patches) = zero    !< Specified x-velocity on patch [m/s].
      real(rk) :: patch_v(max_patches) = zero    !< Specified y-velocity on patch [m/s].
      real(rk) :: patch_w(max_patches) = zero    !< Specified z-velocity on patch [m/s].
      real(rk) :: patch_p(max_patches) = zero    !< Specified static pressure on patch [Pa].
      real(rk) :: patch_dpdn(max_patches) = zero !< Specified pressure gradient on patch [Pa/m].
      real(rk) :: patch_T(max_patches) = 300.0_rk !< Specified temperature on patch [K].

      !> Species Boundary Conditions
      character(len=name_len) :: patch_species_type(max_patches) = "" !< Override BC type for mass fractions.
      real(rk) :: patch_Y(max_species, max_patches) = zero           !< Specified \(Y_k\) mass fractions on patch.

      !> @name Data Output
      character(len=path_len) :: output_dir = "output" !< Directory for result storage.
      logical :: write_vtu = .true.                    !< If true, generates Unstructured VTK files.
      logical :: write_diagnostics = .true.            !< If true, writes global residuals to `diagnostics.csv`.

      !> @name Multi-Species Transport
      logical :: enable_species = .false.                         !< Enable advection-diffusion of mass fractions.
      logical :: enable_reactions = .false.                       !< Reserved for chemical source terms; reactions are not implemented yet.
      integer :: nspecies = 0                                     !< Total number of transport species.
      character(len=name_len) :: species_name(max_species) = ""   !< List of species names.
      real(rk) :: species_diffusivity(max_species) = 0.0_rk       !< Constant species diffusivity \(D_k\) [m^2/s].
      real(rk) :: initial_Y(max_species) = 0.0_rk                 !< Global initial mass fractions in the domain.
      
      !> @name Internal Registry (Reacting discovery)
      integer :: namelist_nspecies = 0                            !< Number of species specified in `case.nml`.
      character(len=name_len) :: namelist_species_name(max_species) = "" !< Names specified in `case.nml`.
      
      !> @name Cantera Bridge Integration
      logical :: enable_cantera_fluid = .false.                   !< Use Cantera for mixture-averaged viscosity; flow density remains `rho`.
      logical :: enable_cantera_species = .false.                 !< Use Cantera for mixture-averaged \(D_k\).
      character(len=path_len) :: cantera_mech_file = "gri30.yaml" !< Path to YAML/CTI mechanism file.
      real(rk) :: background_temp = 300.0_rk                      !< Fixed temperature for property evaluation [K].
      real(rk) :: background_press = 101325.0_rk                  !< Fixed pressure for property evaluation [Pa].

      !> @name Enthalpy Energy Equation Controls
      logical :: enable_energy = .false.                         !< Enable enthalpy/temperature field storage.
      logical :: enable_cantera_thermo = .false.                 !< Use Cantera for sensible h(T,Y,p0) and T(h,Y,p0).
      integer :: thermo_update_interval = 1                     !< Reserved Cantera thermo interval [steps]; must remain 1.
      character(len=name_len) :: thermo_default_species = 'N2'   !< Single-species Cantera thermo fallback when species transport is off.
      real(rk) :: initial_T = 300.0_rk                           !< Initial gas temperature [K].
      real(rk) :: energy_reference_T = 298.15_rk                 !< Reference temperature for sensible enthalpy [K].
      real(rk) :: energy_reference_h = zero                      !< Reference sensible enthalpy for constant-cp mode [J/kg].
      real(rk) :: energy_cp = 1005.0_rk                          !< Constant heat capacity for non-Cantera thermo [J/kg/K].
      real(rk) :: energy_lambda = 2.6e-2_rk                      !< Constant thermal conductivity for non-Cantera thermo [W/m/K].

      !> @name Profiling Controls
      logical :: enable_profiling = .true.                       !< Enable wall-clock profiling.
      logical :: nested_profiling = .true.                       !< If true, print nested profiling tree.
   end type case_params_t

   public :: read_case_params

contains

   !> Orchestrates the reading of all namelist blocks from the configuration file.
   !!
   !! Performs a sequential read of all expected namelist groups and 
   !! triggers a final validation pass to ensure physical consistency.
   !!
   !! @param filename Path to the `.nml` file (usually `case.nml`).
   !! @param params Container to be populated with parsed values.
   subroutine read_case_params(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      call read_mesh_input(filename, params)
      call read_time_input(filename, params)
      call read_fluid_input(filename, params)
      call read_solver_input(filename, params)
      call read_boundary_input(filename, params)
      call read_species_input(filename, params)
      call read_energy_input(filename, params)
      call read_output_input(filename, params)
      call read_profiling_input(filename, params)
      call validate_params(params)
   end subroutine read_case_params


   !> Reads the `&mesh_input` namelist block.
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


   !> Reads the `&time_input` namelist block.
   subroutine read_time_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      integer :: nsteps, output_interval
      real(rk) :: dt, max_cfl
      logical :: use_dynamic_dt
      integer :: unit_id, ios

      namelist /time_input/ nsteps, dt, output_interval, use_dynamic_dt, max_cfl

      nsteps = params%nsteps
      dt = params%dt
      output_interval = params%output_interval
      use_dynamic_dt = params%use_dynamic_dt
      max_cfl = params%max_cfl

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=time_input, iostat=ios)
         close(unit_id)
      end if

      if (ios > 0) call fatal_error('input', 'failed reading &time_input')

      params%nsteps = nsteps
      params%dt = dt
      params%output_interval = output_interval
      params%use_dynamic_dt = use_dynamic_dt
      params%max_cfl = max_cfl
   end subroutine read_time_input


   !> Reads the `&fluid_input` namelist block.
   subroutine read_fluid_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      real(rk) :: rho, nu
      logical :: enable_cantera
      character(len=path_len) :: cantera_mech_file
      real(rk) :: background_temp
      real(rk) :: background_press
      integer :: transport_update_interval
      integer :: unit_id, ios
      logical :: enable_variable_density
      logical :: enable_variable_nu
      namelist /fluid_input/ rho, nu, enable_cantera, cantera_mech_file, background_temp, background_press, transport_update_interval, enable_variable_density, enable_variable_nu

      rho = params%rho
      nu = params%nu
      enable_cantera = params%enable_cantera_fluid
      cantera_mech_file = params%cantera_mech_file
      background_temp = params%background_temp
      background_press = params%background_press
      transport_update_interval = params%transport_update_interval

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
      enable_variable_density = params%enable_variable_density
      enable_variable_nu = params%enable_variable_nu
         read(unit_id, nml=fluid_input, iostat=ios)
         close(unit_id)
      end if

      if (ios > 0) call fatal_error('input', 'failed reading &fluid_input')

      params%rho = rho
      params%nu = nu
      params%enable_cantera_fluid = enable_cantera
         params%enable_variable_density = enable_variable_density
         params%enable_variable_nu = enable_variable_nu
      params%cantera_mech_file = cantera_mech_file
      params%background_temp = background_temp
      params%background_press = background_press
      params%transport_update_interval = transport_update_interval
   end subroutine read_fluid_input


   !> Reads the `&solver_input` namelist block.
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


   !> Reads the `&boundary_input` namelist block.
   subroutine read_boundary_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      integer :: n_patches
      character(len=name_len) :: patch_name(max_patches)
      character(len=name_len) :: patch_type(max_patches)
      character(len=name_len) :: patch_velocity_type(max_patches)
      character(len=name_len) :: patch_pressure_type(max_patches)
      character(len=name_len) :: patch_temperature_type(max_patches)
      character(len=name_len) :: patch_species_type(max_patches)
      real(rk) :: patch_u(max_patches)
      real(rk) :: patch_v(max_patches)
      real(rk) :: patch_w(max_patches)
      real(rk) :: patch_p(max_patches)
      real(rk) :: patch_dpdn(max_patches)
      real(rk) :: patch_T(max_patches)
      real(rk), save :: patch_Y(max_species, max_patches) = zero
      integer :: unit_id, ios, i

      namelist /boundary_input/ n_patches, patch_name, patch_type, &
                                 patch_velocity_type, patch_pressure_type, &
                                 patch_temperature_type, patch_species_type, &
                                 patch_u, patch_v, patch_w, patch_p, patch_dpdn, patch_T, patch_Y

      n_patches = params%n_patches
      patch_name = params%patch_name
      patch_type = params%patch_type
      patch_velocity_type = params%patch_velocity_type
      patch_pressure_type = params%patch_pressure_type
      patch_temperature_type = params%patch_temperature_type
      patch_species_type = params%patch_species_type
      patch_u = params%patch_u
      patch_v = params%patch_v
      patch_w = params%patch_w
      patch_p = params%patch_p
      patch_dpdn = params%patch_dpdn
      patch_T = params%patch_T
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
      params%patch_temperature_type = patch_temperature_type
      params%patch_species_type = patch_species_type

      do i = 1, max_patches
         params%patch_type(i) = trim(lowercase(params%patch_type(i)))
         params%patch_velocity_type(i) = trim(lowercase(params%patch_velocity_type(i)))
         params%patch_pressure_type(i) = trim(lowercase(params%patch_pressure_type(i)))
         params%patch_temperature_type(i) = trim(lowercase(params%patch_temperature_type(i)))
         params%patch_species_type(i) = trim(lowercase(params%patch_species_type(i)))
      end do

      params%patch_u = patch_u
      params%patch_v = patch_v
      params%patch_w = patch_w
      params%patch_p = patch_p
      params%patch_dpdn = patch_dpdn
      params%patch_T = patch_T
      params%patch_Y = patch_Y
   end subroutine read_boundary_input


   !> Reads the `&profiling_input` block.
   !!
   !! This parser is intentionally explicit because profiling is optional and
   !! should be robust to block ordering in case.nml.
   subroutine read_profiling_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      character(len=512) :: line
      character(len=256) :: key, value
      integer :: unit_id, ios, eqpos, comment_pos, comma_pos
      logical :: in_block, found_block

      in_block = .false.
      found_block = .false.

      open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) return

      do
         read(unit_id, '(a)', iostat=ios) line
         if (ios /= 0) exit

         comment_pos = index(line, '!')
         if (comment_pos > 0) line = line(:comment_pos-1)

         line = adjustl(line)
         if (len_trim(line) == 0) cycle

         if (.not. in_block) then
            if (index(lowercase(line), '&profiling_input') == 1) then
               in_block = .true.
               found_block = .true.
            end if
            cycle
         end if

         if (index(line, '/') > 0) exit

         eqpos = index(line, '=')
         if (eqpos <= 0) cycle

         key = trim(adjustl(lowercase(line(:eqpos-1))))
         value = trim(adjustl(lowercase(line(eqpos+1:))))

         comma_pos = index(value, ',')
         if (comma_pos > 0) value = value(:comma_pos-1)
         value = trim(value)

         select case (trim(key))
         case ('enable_profiling')
            call parse_logical_value(value, params%enable_profiling, 'enable_profiling')
         case ('nested_profiling')
            call parse_logical_value(value, params%nested_profiling, 'nested_profiling')
         case default
            call fatal_error('input', 'unknown variable in &profiling_input: '//trim(key))
         end select
      end do

      close(unit_id)

      if (found_block .and. in_block .and. ios /= 0) then
         call fatal_error('input', 'unterminated &profiling_input block')
      end if

   contains

      !> Parses a string representation of a logical value into a Fortran logical.
      !!
      !! Supports common formats like '.true.', 'true', 't', and their false counterparts.
      !!
      !! @param value_text The string to parse.
      !! @param output_value The resulting logical value.
      !! @param field_name The name of the field being parsed (used for error reporting).
      subroutine parse_logical_value(value_text, output_value, field_name)
         character(len=*), intent(in) :: value_text
         logical, intent(out) :: output_value
         character(len=*), intent(in) :: field_name

         select case (trim(value_text))
         case ('.true.', 'true', 't', '.t.')
            output_value = .true.
         case ('.false.', 'false', 'f', '.f.')
            output_value = .false.
         case default
            call fatal_error('input', 'invalid logical for '//trim(field_name)//': '//trim(value_text))
         end select
      end subroutine parse_logical_value

   end subroutine read_profiling_input



   !> Reads the `&output_input` namelist block.
   !> Reads the `&energy_input` namelist block.
   !!
   !! Reads controls for sensible-enthalpy transport. Cantera thermo, when
   !! enabled, preserves transported `h` and recovers `T(h,Y,p0)` every energy
   !! step; `thermo_update_interval` is reserved and must remain 1.
   subroutine read_energy_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      logical :: enable_energy
      logical :: enable_cantera_thermo
      integer :: thermo_update_interval
      character(len=name_len) :: thermo_default_species
      real(rk) :: initial_T
      real(rk) :: energy_reference_T
      real(rk) :: energy_reference_h
      real(rk) :: energy_cp
      real(rk) :: energy_lambda
      integer :: unit_id, ios

      namelist /energy_input/ enable_energy, enable_cantera_thermo, thermo_default_species, thermo_update_interval, &
                               initial_T, energy_reference_T, energy_reference_h, &
                               energy_cp, energy_lambda

      enable_energy = params%enable_energy
      enable_cantera_thermo = params%enable_cantera_thermo
      thermo_update_interval = params%thermo_update_interval
      thermo_default_species = params%thermo_default_species
      initial_T = params%initial_T
      energy_reference_T = params%energy_reference_T
      energy_reference_h = params%energy_reference_h
      energy_cp = params%energy_cp
      energy_lambda = params%energy_lambda

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=energy_input, iostat=ios)
         close(unit_id)
      end if

      if (ios /= 0 .and. ios /= -1) then
         call fatal_error('input', 'failed reading &energy_input. Check for unknown variables or typos.')
      end if

      if (ios == 0) then
         params%enable_energy = enable_energy
         params%enable_cantera_thermo = enable_cantera_thermo
         params%thermo_update_interval = thermo_update_interval
         params%thermo_default_species = trim(thermo_default_species)
         params%initial_T = initial_T
         params%energy_reference_T = energy_reference_T
         params%energy_reference_h = energy_reference_h
         params%energy_cp = energy_cp
         params%energy_lambda = energy_lambda
      end if
   end subroutine read_energy_input


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


   !> Reads the `&species_input` namelist block.
   !!
   !! Implements strict error checking to catch typos in chemical 
   !! species names or advection settings.
   subroutine read_species_input(filename, params)
      character(len=*), intent(in) :: filename
      type(case_params_t), intent(inout) :: params

      logical :: enable_species, enable_reactions, enable_cantera
      integer :: nspecies
      character(len=name_len), save :: species_name(max_species) = ""
      real(rk), save :: species_diffusivity(max_species) = zero
      real(rk), save :: initial_Y(max_species) = zero
      integer :: unit_id, ios

      namelist /species_input/ enable_species, enable_reactions, enable_cantera, &
                                nspecies, species_name, species_diffusivity, initial_Y

      enable_species = params%enable_species
      enable_reactions = params%enable_reactions
      enable_cantera = params%enable_cantera_species
      nspecies = params%nspecies
      species_name = params%species_name
      species_diffusivity = params%species_diffusivity
      initial_Y = params%initial_Y

      call open_namelist_file(filename, unit_id, ios)

      if (ios == 0) then
         read(unit_id, nml=species_input, iostat=ios)
         close(unit_id)
      end if

      if (ios /= 0 .and. ios /= -1) then
         call fatal_error('input', 'failed reading &species_input. Check for unknown variables or typos.')
      end if

      if (ios == 0) then
         params%enable_species = enable_species
         params%enable_reactions = enable_reactions
         params%enable_cantera_species = enable_cantera
         params%nspecies = nspecies
         params%species_name = species_name
         params%species_diffusivity = species_diffusivity
         params%initial_Y = initial_Y
         
         ! Save original namelist values
         params%namelist_nspecies = nspecies
         params%namelist_species_name = species_name
      end if
   end subroutine read_species_input


   !> Helper routine to safely open a namelist file for reading.
   subroutine open_namelist_file(filename, unit_id, ios)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: unit_id
      integer, intent(out) :: ios

      open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)

      if (ios /= 0) unit_id = -1
   end subroutine open_namelist_file


   !> Validates all parsed parameters against physical and algorithmic limits.
   !!
   !! Triggers `fatal_error` for non-physical values (e.g., negative density, 
   !! zero timestep) or unsupported configurations.
   subroutine validate_params(params)
      type(case_params_t), intent(in) :: params

      if (params%nsteps < 0) call fatal_error('input', 'nsteps must be non-negative')
      if (params%dt <= zero) call fatal_error('input', 'dt must be positive')
      if (params%output_interval <= 0) call fatal_error('input', 'output_interval must be positive')
      if (params%rho <= zero) call fatal_error('input', 'rho must be positive')
      if (params%enable_variable_density) then
         call fatal_error('input', 'enable_variable_density=.true. requested, but variable-density flow is not implemented yet; rho_thermo is diagnostic only')
      end if
      if (params%nu < zero) call fatal_error('input', 'nu must be non-negative')
      if (params%enable_variable_nu .and. .not. params%enable_cantera_fluid) then
         call fatal_error('input', 'enable_variable_nu=.true. requires fluid_input enable_cantera=.true.; otherwise use constant nu')
      end if
      if (params%transport_update_interval <= 0) call fatal_error('input', 'transport_update_interval must be positive')
      if (params%pressure_max_iter <= 0) call fatal_error('input', 'pressure_max_iter must be positive')
      if (params%pressure_tol <= zero) call fatal_error('input', 'pressure_tol must be positive')

      if (params%initial_T <= zero) call fatal_error('input', 'initial_T must be positive')
      if (params%energy_reference_T <= zero) call fatal_error('input', 'energy_reference_T must be positive')
      if (params%energy_cp <= zero) call fatal_error('input', 'energy_cp must be positive')
      if (params%energy_lambda < zero) call fatal_error('input', 'energy_lambda must be non-negative')
      if (params%thermo_update_interval <= 0) call fatal_error('input', 'thermo_update_interval must be positive')
      if (params%enable_cantera_thermo .and. params%thermo_update_interval /= 1) then
         call fatal_error('input', 'thermo_update_interval values other than 1 are reserved; Cantera thermo is currently updated every energy step')
      end if
      if (params%enable_cantera_thermo .and. len_trim(params%thermo_default_species) == 0) &
         call fatal_error('input', 'thermo_default_species must not be empty when enable_cantera_thermo is true')

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


   !> Ensures all boundary patch names and types are non-empty.
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

         if (params%enable_energy .and. params%patch_T(i) <= zero) then
            call fatal_error('input', 'patch_T entry must be positive when energy is enabled')
         end if
      end do
   end subroutine validate_boundary_arrays

end module mod_input