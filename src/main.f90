!> Main entry point for the LowMachReact-Hex constant-density FV solver.
!!
!! This program orchestrates the entire simulation lifecycle:
!! 1. **Initialization**: Starts MPI, parses the `case.nml` namelist, and reads the mesh.
!! 2. **Domain Decomposition**: Sets up the replicated mesh MPI ranks for flow and radiation.
!! 3. **Field Setup**: Allocates and initializes velocity, pressure, species, and energy fields.
!! 4. **Time Integration Loop**: Executes projection, species transport, and sensible-enthalpy transport.
!! 5. **Diagnostics & Output**: Computes global observables and writes VTU/CSV files.
!! 6. **Finalization**: Safely releases all allocated memory and shuts down MPI.
program lowmach_react_hex
   use mpi_f08
   use mod_kinds, only : rk, zero, one, output_unit, fatal_error
   use mod_input, only : case_params_t, read_case_params
   use mod_mesh_types, only : mesh_t, mesh_finalize
   use mod_mesh_io, only : read_native_mesh
   use mod_mpi_flow, only : flow_mpi_t, mpi_flow_startup, mpi_flow_shutdown, &
                            flow_mpi_initialize, flow_mpi_finalize
   use mod_mpi_radiation, only : radiation_mpi_t, radiation_mpi_initialize, &
                                 radiation_mpi_finalize
   use mod_bc, only : bc_set_t, build_bc_set, finalize_bc_set
   use mod_fields, only : flow_fields_t, initialize_fields, finalize_fields
   use mod_flow_projection, only : solver_stats_t, advance_projection_step, compute_flow_diagnostics, compute_and_update_cfl
   use mod_output, only : prepare_output, write_diagnostics_header, write_diagnostics_row, &
                          write_vtu_unstructured, write_mesh_summary, write_pvd_collection
   use mod_species, only : species_fields_t, initialize_species, finalize_species, &
                           advance_species_transport
   use mod_energy, only : energy_fields_t, initialize_energy, finalize_energy, &
                           advance_energy_transport, &
                           write_energy_diagnostics_header, write_energy_diagnostics_row
   use mod_transport_properties, only : transport_properties_t, initialize_transport, &
                                        finalize_transport, update_transport_properties
   use mod_profiler, only : profiler_start, profiler_stop, profiler_report, profiler_configure
   implicit none

   type(case_params_t) :: params    !< Parsed simulation parameters.
   type(mesh_t) :: mesh             !< Global computational mesh.
   type(flow_mpi_t) :: flow_mpi     !< MPI context for flow solver.
   type(radiation_mpi_t) :: rad_mpi !< MPI context for radiation solver.
   type(bc_set_t) :: bc             !< Boundary condition definitions.
   type(flow_fields_t) :: fields    !< Hydrodynamic fields (U, P).
   type(solver_stats_t) :: stats    !< Runtime diagnostics and solver performance metrics.
   type(species_fields_t) :: species !< Species mass fractions.
   type(energy_fields_t) :: energy   !< Enthalpy/temperature/radiation-source fields.
   type(transport_properties_t) :: transport !< Physical properties (rho, mu, Dk).

   character(len=256) :: case_file
   integer :: step, ierr, world_rank
   logical :: do_output_step
   real(rk) :: time, t0

   ! 1. Initialize MPI and read runtime configuration.
   call mpi_flow_startup()
   call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
   if (ierr /= MPI_SUCCESS) call fatal_error('main', 'MPI_Comm_rank failed')

   t0 = real(mpi_wtime(), rk)

   call get_case_filename(case_file)
   call read_case_params(trim(case_file), params)
   call profiler_configure(params%enable_profiling, params%nested_profiling)

   if (world_rank == 0) then
      write(output_unit,'(a,l1,a,l1)') 'profiling config: enable=', &
         params%enable_profiling, ' nested=', params%nested_profiling
   end if

   call profiler_start('Total_Simulation')

   ! 2. Setup simulation environment and read input data.
   call read_native_mesh(trim(params%mesh_dir), mesh)
   
   ! 3. Initialize parallel contexts.
   call initialize_transport(mesh, params, transport)
   call flow_mpi_initialize(mesh, flow_mpi, MPI_COMM_WORLD, max(4, params%nspecies))
   call radiation_mpi_initialize(rad_mpi, MPI_COMM_WORLD)
   
   ! 4. Prepare physical fields and BCs.
   call build_bc_set(mesh, params, bc)
   call initialize_fields(mesh, fields)
   if (params%enable_species) then
      call initialize_species(mesh, params, species)
   end if
   if (params%enable_energy) then
      if (params%enable_species) then
         call initialize_energy(mesh, params, energy, species%Y)
      else
         call initialize_energy(mesh, params, energy)
      end if
   end if

   if (flow_mpi%rank == 0 .and. params%enable_species) then
      print *, "species: ", species%nspecies
   end if

   ! 5. Setup output files.
   call prepare_output(params, flow_mpi)
   call write_mesh_summary(params, flow_mpi, mesh)
   call write_diagnostics_header(params, flow_mpi)
   call write_energy_diagnostics_header(params, flow_mpi)
   call write_pvd_collection(params, flow_mpi, params%nsteps, params%output_interval, params%dt)

   ! 6. Initial diagnostics and visualization snapshot.
   time = zero
   call compute_flow_diagnostics(mesh, flow_mpi, bc, params, fields, stats)
   call write_diagnostics_row(params, flow_mpi, 0, time, stats)
   call write_energy_diagnostics_row(mesh, flow_mpi, params, energy, 0, time)
   call write_vtu_unstructured(params, flow_mpi, mesh, fields, species, energy, transport, 0)

   if (flow_mpi%rank == 0) then
      write(output_unit,'(a)') 'LowMachReact-Hex hexahedral low-Mach FV solver'
      write(output_unit,'(a,a)') 'case: ', trim(case_file)
      write(output_unit,'(a,i0)') 'cells: ', mesh%ncells
      write(output_unit,'(a,i0)') 'flow/chemistry MPI ranks: ', flow_mpi%nprocs
      write(output_unit,'(a,i0)') 'radiation MPI ranks: ', rad_mpi%nprocs
      write(output_unit,'(a,es12.5,a)') 'Flow density mode: constant rho = ', params%rho, ' kg/m^3'
      if (params%enable_variable_nu) then
         write(output_unit,'(a)') 'Flow viscosity mode: variable Cantera mu/nu enabled'
      else
         write(output_unit,'(a,es12.5,a)') 'Flow viscosity mode: constant nu = ', params%nu, ' m^2/s'
      end if
      if (params%enable_energy .and. params%enable_cantera_thermo) then
         write(output_unit,'(a)') 'Cantera rho_thermo: diagnostic only, not used by projection'
      end if
      if (params%enable_cantera_fluid .or. params%enable_cantera_species) then
         write(output_unit,'(a,i0,a)') 'Cantera transport update interval: ', &
            params%transport_update_interval, ' step(s) [mu/D_k only]'
         if (params%enable_energy) then
            write(output_unit,'(a)') 'Cantera transport temperature source: energy%T'
         else
            write(output_unit,'(a,es12.5,a)') 'Cantera transport temperature source: background_temp = ', &
               params%background_temp, ' K'
         end if
      end if
      if (params%enable_energy .and. params%enable_cantera_thermo) then
         write(output_unit,'(a)') 'Cantera thermo update interval: every energy step'
         write(output_unit,'(a,es12.5,a)') 'Cantera thermodynamic pressure p0: ', &
            params%background_press, ' Pa'
      end if
   end if

   ! 7. Main Time-Stepping Loop
   do step = 1, params%nsteps
      do_output_step = (mod(step, params%output_interval) == 0) .or. step == params%nsteps

      ! A. Adaptive time-stepping / CFL diagnostics.
      ! Dynamic dt requires CFL every step. For fixed-dt runs, compute CFL only
      ! when diagnostics are written.
      if (params%use_dynamic_dt .or. do_output_step) then
         call profiler_start('CFL_Update')
         call compute_and_update_cfl(mesh, flow_mpi, params, fields, stats)
         call profiler_stop('CFL_Update')
      end if

      ! B. Dynamic transport-property update via Cantera.
      ! transport_update_interval gates Cantera mu/D_k refresh; enable_variable_nu controls whether mu affects flow; flow density remains params%rho.
      ! Cantera h/T thermo recovery remains inside the energy update path.
      if (mod(step-1, params%transport_update_interval) == 0 .or. step == 1) then
         call profiler_start('Transport_Update')
         if (params%enable_species) then
            if (params%enable_energy) then
               call update_transport_properties(mesh, flow_mpi, params, species%Y, transport, T_state=energy%T)
            else
               call update_transport_properties(mesh, flow_mpi, params, species%Y, transport)
            end if
         else
            if (params%enable_energy) then
               call update_transport_properties(mesh, flow_mpi, params, transport=transport, T_state=energy%T)
            else
               call update_transport_properties(mesh, flow_mpi, params, transport=transport)
            end if
         end if
         call profiler_stop('Transport_Update')
      end if
      
      ! C. Advance Momentum & Pressure (Projection Method).
      call profiler_start('Projection_Step')
      call advance_projection_step(mesh, flow_mpi, bc, params, transport, fields, stats)
      call profiler_stop('Projection_Step')
      
      ! D. Advance Species Transport.
      if (params%enable_species) then
         call profiler_start('Species_Transport')
         call advance_species_transport(mesh, flow_mpi, bc, params, fields, species, transport)
         call profiler_stop('Species_Transport')
      end if
      
      ! E. Advance passive sensible-enthalpy transport.
      if (params%enable_energy) then
         call profiler_start('Energy_Transport')
         if (params%enable_species) then
            call advance_energy_transport(mesh, flow_mpi, bc, params, fields, energy, species%Y)
         else
            call advance_energy_transport(mesh, flow_mpi, bc, params, fields, energy)
         end if
         call profiler_stop('Energy_Transport')
      end if
      
      time = time + params%dt

      ! E. Periodic Output and Diagnostics.
      if (do_output_step) then
         stats%wall_time = real(mpi_wtime(), rk) - t0
         call profiler_start('Flow_Diagnostics')
         call compute_global_observables(mesh, flow_mpi, fields, species, transport, stats)
         call profiler_stop('Flow_Diagnostics')

         call profiler_start('Diagnostics_Write_Flow')
         call write_diagnostics_row(params, flow_mpi, step, time, stats)
         call profiler_stop('Diagnostics_Write_Flow')

         if (params%enable_energy) then
            call profiler_start('Diagnostics_Write_Energy')
            call write_energy_diagnostics_row(mesh, flow_mpi, params, energy, step, time)
            call profiler_stop('Diagnostics_Write_Energy')
         end if

         call profiler_start('Output_Write_VTU')
         call write_vtu_unstructured(params, flow_mpi, mesh, fields, species, energy, transport, step)
         call profiler_stop('Output_Write_VTU')

         if (flow_mpi%rank == 0) then
            write(output_unit,'(a,i0,2x,a,es12.5,2x,a,f8.2,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5,2x,a,i0,2x,a,es12.5)') &
               'step=', step, 'time=', time, 'clock=', stats%wall_time, 'dt=', params%dt, 'max_div=', stats%max_divergence, &
               'cfl=', stats%cfl, 'ke=', stats%kinetic_energy, 'piter=', stats%pressure_iterations, '|U|max=', stats%max_velocity
         end if
      end if
   end do

   ! 8. Clean up and finalize simulation.
   call finalize_bc_set(bc)
   call finalize_fields(fields)
   if (params%enable_species) then
      call finalize_species(species)
   end if
   if (params%enable_energy) then
      call finalize_energy(energy)
   end if
   call finalize_transport(transport)
   call mesh_finalize(mesh)

   call profiler_stop('Total_Simulation')
   call profiler_report(flow_mpi%comm, flow_mpi%rank, flow_mpi%nprocs)

   call radiation_mpi_finalize(rad_mpi)
   call flow_mpi_finalize(flow_mpi)
   call mpi_flow_shutdown()

contains

   !> Computes global integral quantities for diagnostic reporting.
   !!
   !! Performs Allreduce operations to find the maximum velocity magnitude, 
   !! total system mass, and minimum species mass fraction.
   subroutine compute_global_observables(mesh, flow, fields, species, transport, stats)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(flow_fields_t), intent(in) :: fields
      type(species_fields_t), intent(in) :: species
      type(transport_properties_t), intent(in) :: transport
      type(solver_stats_t), intent(inout) :: stats

      integer :: c
      real(rk) :: local_max_vel, local_min_y, local_total_mass
      real(rk) :: global_max_vel, global_min_y, global_total_mass
      real(rk) :: u, v, w, vel_mag

      local_max_vel = zero
      local_min_y = huge(one)
      local_total_mass = zero

      do c = 1, mesh%ncells
         if (.not. flow%owned(c)) cycle

         u = fields%u(1, c)
         v = fields%u(2, c)
         w = fields%u(3, c)
         vel_mag = sqrt(u**2 + v**2 + w**2)
         if (vel_mag > local_max_vel) local_max_vel = vel_mag

         local_total_mass = local_total_mass + (transport%rho(c) * mesh%cells(c)%volume)

         if (species%nspecies > 0) then
            local_min_y = min(local_min_y, minval(species%Y(:, c)))
         end if
      end do

      call MPI_Allreduce(local_max_vel, global_max_vel, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm)
      call MPI_Allreduce(local_total_mass, global_total_mass, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm)
      
      if (species%nspecies > 0) then
         call MPI_Allreduce(local_min_y, global_min_y, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm)
      else
         global_min_y = zero
      end if

      stats%max_velocity = global_max_vel
      stats%total_mass = global_total_mass
      stats%min_species_y = global_min_y

   end subroutine compute_global_observables


   !> Parses command line arguments to find the case configuration file.
   !!
   !! If no argument is provided, defaults to `case.nml`.
   subroutine get_case_filename(filename)
      character(len=*), intent(out) :: filename

      integer :: argc

      argc = command_argument_count()
      if (argc >= 1) then
         call get_command_argument(1, filename)
      else
         filename = 'case.nml'
      end if

      if (len_trim(filename) == 0) call fatal_error('main', 'empty case filename')
   end subroutine get_case_filename

end program lowmach_react_hex
