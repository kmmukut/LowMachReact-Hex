!> Main entry point for the LowMachReact-Hex hexahedral low-Mach FV solver.
program lowmach_react_hex
   use mpi_f08
   use mod_kinds, only : rk, zero, output_unit, fatal_error
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
   use mod_transport_properties, only : transport_properties_t, initialize_transport, &
                                        finalize_transport, update_transport_properties
   implicit none

   type(case_params_t) :: params
   type(mesh_t) :: mesh
   type(flow_mpi_t) :: flow_mpi
   type(radiation_mpi_t) :: rad_mpi
   type(bc_set_t) :: bc
   type(flow_fields_t) :: fields
   type(solver_stats_t) :: stats
   type(species_fields_t) :: species
   type(transport_properties_t) :: transport

   character(len=256) :: case_file
   integer :: step
   real(rk) :: time

   call mpi_flow_startup()

   call get_case_filename(case_file)
   call read_case_params(trim(case_file), params)
   call read_native_mesh(trim(params%mesh_dir), mesh)
   call flow_mpi_initialize(mesh, flow_mpi, MPI_COMM_WORLD)
   call radiation_mpi_initialize(rad_mpi, MPI_COMM_WORLD)
   call build_bc_set(mesh, params, bc)
   call initialize_fields(mesh, fields)
   call initialize_transport(mesh, params, transport)
   if (params%enable_species) then
      call initialize_species(mesh, params, species)
   end if

   if (flow_mpi%rank == 0 .and. params%enable_species) then
      print *, "species: ", species%nspecies
   end if

   call prepare_output(params, flow_mpi)
   call write_mesh_summary(params, flow_mpi, mesh)
   call write_diagnostics_header(params, flow_mpi)
   call write_pvd_collection(params, flow_mpi, params%nsteps, params%output_interval, params%dt)

   time = zero
   call compute_flow_diagnostics(mesh, flow_mpi, bc, params, fields, stats)
   call write_diagnostics_row(params, flow_mpi, 0, time, stats)
   call write_vtu_unstructured(params, flow_mpi, mesh, fields, species, 0)

   if (flow_mpi%rank == 0) then
      write(output_unit,'(a)') 'LowMachReact-Hex hexahedral low-Mach FV solver'
      write(output_unit,'(a,a)') 'case: ', trim(case_file)
      write(output_unit,'(a,i0)') 'cells: ', mesh%ncells
      write(output_unit,'(a,i0)') 'flow/chemistry MPI ranks: ', flow_mpi%nprocs
      write(output_unit,'(a,i0)') 'radiation MPI ranks: ', rad_mpi%nprocs
   end if

   do step = 1, params%nsteps
      call compute_and_update_cfl(mesh, flow_mpi, params, fields, stats)

      if (mod(step-1, params%transport_update_interval) == 0 .or. step == 1) then
         if (params%enable_species) then
            call update_transport_properties(mesh, params, species%Y, transport)
         else
            call update_transport_properties(mesh, params, transport=transport)
         end if
      end if
      call advance_projection_step(mesh, flow_mpi, bc, params, transport, fields, stats)
      if (params%enable_species) then
         call advance_species_transport(mesh, flow_mpi, bc, params, fields, species, transport)
      end if
      time = time + params%dt

      if (mod(step, params%output_interval) == 0 .or. step == params%nsteps) then
         call write_diagnostics_row(params, flow_mpi, step, time, stats)
         call write_vtu_unstructured(params, flow_mpi, mesh, fields, species, step)

         if (flow_mpi%rank == 0) then
            write(output_unit,'(a,i0,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5,2x,a,es12.5,2x,a,i0)') &
               'step=', step, 'time=', time, 'dt=', params%dt, 'max_div=', stats%max_divergence, &
               'cfl=', stats%cfl, 'ke=', stats%kinetic_energy, 'piter=', stats%pressure_iterations
         end if
      end if
   end do

   call finalize_transport(transport)
   if (params%enable_species) then
      call finalize_species(species)
   end if
   call finalize_fields(fields)
   call finalize_bc_set(bc)
   call radiation_mpi_finalize(rad_mpi)
   call flow_mpi_finalize(flow_mpi)
   call mesh_finalize(mesh)
   call mpi_flow_shutdown()

contains

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
