!> Output management for VTK visualization and solver diagnostics.
!!
!! This module handles the generation of simulation results in modern XML-based 
!! VTK format (.vtu) and CSV-based global diagnostics. It manages the 
!! creation of the output directory, writing mesh summaries, and 
!! generating PVD collection files for time-series visualization in ParaView.
!! Most arrays are globally allocated, but VTU piece files write owned cells
!! only. Variable-density projection validation should use
!! `divu_minus_S_projection` diagnostics, not raw `divergence`.
module mod_output
   use mpi_f08
   use mod_kinds, only : rk, zero, path_len, fatal_error
   use mod_input, only : case_params_t
   use mod_mesh_types, only : mesh_t
   use mod_mpi_flow, only : flow_mpi_t, flow_gather_owned_scalar_root, &
                            flow_gather_owned_matrix_root
   use mod_fields, only : flow_fields_t
   use mod_flow_projection, only : solver_stats_t
   use mod_species, only : species_fields_t
   use mod_energy, only : energy_fields_t
   use mod_transport_properties, only : transport_properties_t
   implicit none

   private

   public :: prepare_output, write_diagnostics_header, write_diagnostics_row
   public :: write_vtu_unstructured, write_mesh_summary, write_pvd_collection, write_variable_density_diagnostics
   public :: write_species_energy_conservation_diagnostics
   public :: write_enthalpy_energy_budget_diagnostics

contains

   !> Creates the output directory specified in the case parameters.

   subroutine prepare_output(params, flow)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow

      integer :: exitstat
      character(len=path_len + 32) :: command

      if (flow%rank /= 0) return

      command = 'mkdir -p ' // trim(params%output_dir)
      call execute_command_line(trim(command), exitstat=exitstat)
      if (exitstat /= 0) call fatal_error('output', 'failed to create output directory: ' // trim(params%output_dir))

      command = 'mkdir -p ' // trim(params%output_dir) // '/VTK'
      call execute_command_line(trim(command), exitstat=exitstat)
      if (exitstat /= 0) call fatal_error('output', 'failed to create VTK output directory: ' // trim(params%output_dir) // '/VTK')

      command = 'mkdir -p ' // trim(params%output_dir) // '/diagnostics'
      call execute_command_line(trim(command), exitstat=exitstat)
      if (exitstat /= 0) call fatal_error('output', 'failed to create diagnostics output directory: ' // trim(params%output_dir) // '/diagnostics')

   end subroutine prepare_output



   !> Writes the CSV header for global simulation diagnostics.
   subroutine write_diagnostics_header(params, flow)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow

      integer :: unit_id
      character(len=path_len + 32) :: filename

      if (flow%rank /= 0 .or. .not. params%write_diagnostics) return

      filename = trim(params%output_dir)//'/diagnostics/diagnostics.csv'
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a)') 'step,time,dt,max_divergence,rms_divergence,net_boundary_flux,'// &
                           'kinetic_energy,pressure_iterations,pressure_iterations_total,'// &
                           'pressure_iterations_max,pressure_iterations_avg,pressure_solve_count,'// &
                           'pressure_residual,cfl,wall_time,max_velocity,total_mass,min_species_y'

      close(unit_id)
   end subroutine write_diagnostics_header


   !> Appends a new row of diagnostic data to the CSV file.
   subroutine write_diagnostics_row(params, flow, step, time, stats)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow
      integer, intent(in) :: step
      real(rk), intent(in) :: time
      type(solver_stats_t), intent(in) :: stats

      integer :: unit_id
      character(len=path_len + 32) :: filename

      if (flow%rank /= 0 .or. .not. params%write_diagnostics) return

      filename = trim(params%output_dir)//'/diagnostics/diagnostics.csv'
      open(newunit=unit_id, file=trim(filename), status='old', position='append', action='write')

      write(unit_id,'(i0,a,ES26.16E4,a,ES26.16E4,a,ES26.16E4,a,ES26.16E4,a,ES26.16E4,a,ES26.16E4,a,'// &
                    'i0,a,i0,a,i0,a,ES26.16E4,a,i0,a,ES26.16E4,a,ES26.16E4,a,ES26.16E4,a,'// &
                    'ES26.16E4,a,ES26.16E4,a,ES26.16E4)') &
            step, ',', time, ',', params%dt, ',', stats%max_divergence, ',', &
            stats%rms_divergence, ',', stats%net_boundary_flux, ',', &
            stats%kinetic_energy, ',', stats%pressure_iterations, ',', &
            stats%pressure_iterations_total, ',', stats%pressure_iterations_max, ',', &
            stats%pressure_iterations_avg, ',', stats%pressure_solve_count, ',', &
            stats%pressure_residual, ',', stats%cfl, ',', stats%wall_time, ',', &
            stats%max_velocity, ',', stats%total_mass, ',', stats%min_species_y

      close(unit_id)
   end subroutine write_diagnostics_row


   !> Writes a human-readable summary of the mesh connectivity and patches.
   subroutine write_mesh_summary(params, flow, mesh)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow
      type(mesh_t), intent(in) :: mesh

      integer :: unit_id, p
      character(len=path_len + 32) :: filename

      if (flow%rank /= 0) return

      filename = trim(params%output_dir)//'/mesh_summary.txt'
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a,i0)') 'points: ', mesh%npoints
      write(unit_id,'(a,i0)') 'cells: ', mesh%ncells
      write(unit_id,'(a,i0)') 'faces: ', mesh%nfaces
      write(unit_id,'(a,i0)') 'patches: ', mesh%npatches

      do p = 1, mesh%npatches
         write(unit_id,'(i0,1x,a,1x,i0)') mesh%patches(p)%id, trim(mesh%patches(p)%name), mesh%patches(p)%nfaces
      end do

      close(unit_id)
   end subroutine write_mesh_summary


   !> Writes the full flow field to an XML Unstructured Grid file (.vtu).
   !!
   !! Fields included:
   !! - Velocity (Cell Vector)
   !! - Pressure (Cell Scalar)
   !! - Flow density rho (Cell Scalar)
   !! - Kinematic viscosity nu (Cell Scalar)
   !! - Divergence (Cell Scalar)
   !! - Species Mass Fractions (Cell Scalars, if enabled)
   subroutine write_vtu_unstructured(params, flow, mesh, fields, species, energy, transport, step)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(inout) :: flow
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(in) :: fields
      type(species_fields_t), intent(in) :: species
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: step

      integer :: unit_id
      integer :: p, c, n, m, k, n_owned
      character(len=path_len + 64) :: filename, local_filename
      integer, allocatable :: owned_indices(:)

      if (.not. params%write_vtu) return

      ! Count owned cells and collect their indices
      n_owned = 0
      do c = 1, mesh%ncells
         if (flow%owned(c)) n_owned = n_owned + 1
      end do

      allocate(owned_indices(n_owned))
      n_owned = 0
      do c = 1, mesh%ncells
         if (flow%owned(c)) then
            n_owned = n_owned + 1
            owned_indices(n_owned) = c
         end if
      end do

      ! Every rank writes its own piece file
      write(local_filename,'("flow_",i6.6,"_P",i4.4,".vtu")') step, flow%rank
      filename = trim(params%output_dir)//'/VTK/'//trim(local_filename)
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a)') '<?xml version="1.0"?>'
      write(unit_id,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
      write(unit_id,'(a)') '  <UnstructuredGrid>'
      write(unit_id,'(a,i0,a,i0,a)') '    <Piece NumberOfPoints="', mesh%npoints, &
                                      '" NumberOfCells="', n_owned, '">'

      ! ------------------------------------------------------------
      ! Correct XML VTK Piece order:
      !
      ! PointData
      ! CellData
      ! Points
      ! Cells
      ! ------------------------------------------------------------

      write(unit_id,'(a)') '      <PointData>'
      write(unit_id,'(a)') '      </PointData>'

      write(unit_id,'(a)') '      <CellData Scalars="pressure" Vectors="velocity">'

      write(unit_id,'(a)') '        <DataArray type="Float64" Name="velocity" NumberOfComponents="3" format="ascii">'
      do n = 1, n_owned
         c = owned_indices(n)
         write(unit_id,'(3(ES26.16E4,1x))') fields%u(1,c), fields%u(2,c), fields%u(3,c)
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a,a,a)') '        <DataArray type="Float64" Name="pressure" format="ascii">'
      do n = 1, n_owned
         c = owned_indices(n)
         write(unit_id,'(ES26.16E4)') fields%p(c)
      end do
      write(unit_id,'(a)') '        </DataArray>'


      write(unit_id,'(a)') '        <DataArray type="Float64" Name="thermo_pressure" format="ascii">'
      do n = 1, n_owned
         write(unit_id,'(ES26.16E4)') params%background_press
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a,a,a)') '        <DataArray type="Float64" Name="divergence" format="ascii">'
      do n = 1, n_owned
         c = owned_indices(n)
         write(unit_id,'(ES26.16E4)') fields%div(c)
      end do
      write(unit_id,'(a)') '        </DataArray>'

      if (.not. allocated(transport%rho) .or. .not. allocated(transport%nu)) then
         call fatal_error('output', 'transport rho/nu output requested but transport arrays are not allocated')
      end if

      write(unit_id,'(a)') '        <DataArray type="Float64" Name="rho" format="ascii">'
      do n = 1, n_owned
         c = owned_indices(n)
         write(unit_id,'(ES26.16E4)') transport%rho(c)
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a)') '        <DataArray type="Float64" Name="nu" format="ascii">'
      do n = 1, n_owned
         c = owned_indices(n)
         write(unit_id,'(ES26.16E4)') transport%nu(c)
      end do
      write(unit_id,'(a)') '        </DataArray>'

      if (params%enable_energy) then
         if (.not. allocated(energy%T) .or. .not. allocated(energy%h) .or. &
             .not. allocated(energy%qrad)) then
            call fatal_error('output', 'energy output requested but energy arrays are not allocated')
         end if

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="temperature" format="ascii">'
         do n = 1, n_owned
            c = owned_indices(n)
            write(unit_id,'(ES26.16E4)') energy%T(c)
         end do
         write(unit_id,'(a)') '        </DataArray>'

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="enthalpy" format="ascii">'
         do n = 1, n_owned
            c = owned_indices(n)
            write(unit_id,'(ES26.16E4)') energy%h(c)
         end do
         write(unit_id,'(a)') '        </DataArray>'

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="qrad" format="ascii">'
         do n = 1, n_owned
            c = owned_indices(n)
            write(unit_id,'(ES26.16E4)') energy%qrad(c)
         end do
         write(unit_id,'(a)') '        </DataArray>'


         if (params%enable_species_enthalpy_diffusion .and. allocated(energy%species_enthalpy_diffusion)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="species_enthalpy_diffusion" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(ES26.16E4)') energy%species_enthalpy_diffusion(c)
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if

         if (allocated(energy%cp)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="cp" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(ES26.16E4)') energy%cp(c)
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if

         if (allocated(energy%lambda)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="thermal_conductivity" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(ES26.16E4)') energy%lambda(c)
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if


         if (allocated(energy%lambda) .and. allocated(energy%cp)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="thermal_diffusivity" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(ES26.16E4)') energy%lambda(c) / max(params%rho * energy%cp(c), tiny(1.0_rk))
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if

         if (allocated(energy%rho_thermo)) then
            write(unit_id,'(a)') '        <DataArray type="Float64" Name="rho_thermo" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(ES26.16E4)') energy%rho_thermo(c)
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end if

         if (params%enable_variable_density) then
            call write_energy_reconciliation_vtu_arrays(unit_id, mesh, flow, energy, transport)
         end if

      end if

      if (params%enable_species .and. params%nspecies > 0) then

         do k = 1, params%nspecies
            write(unit_id,'(a,a,a)') '        <DataArray type="Float64" Name="Y_', &
               trim(params%species_name(k)), '" format="ascii">'
            do n = 1, n_owned
               c = owned_indices(n)
               write(unit_id,'(ES26.16E4)') species%Y(k,c)
            end do
            write(unit_id,'(a)') '        </DataArray>'
         end do

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="sum_Y" format="ascii">'
         do n = 1, n_owned
            c = owned_indices(n)
            write(unit_id,'(ES26.16E4)') sum(species%Y(:,c))
         end do
         write(unit_id,'(a)') '        </DataArray>'

         if (allocated(transport%diffusivity)) then
            do k = 1, params%nspecies
               write(unit_id,'(a,a,a)') '        <DataArray type="Float64" Name="D_', &
                  trim(params%species_name(k)), '" format="ascii">'
               do n = 1, n_owned
                  c = owned_indices(n)
                  write(unit_id,'(ES26.16E4)') transport%diffusivity(k,c)
               end do
               write(unit_id,'(a)') '        </DataArray>'
            end do
         end if

      end if

      ! Helpful debug scalar so you can color by cell_id in ParaView.
      write(unit_id,'(a)') '        <DataArray type="Int32" Name="cell_id" format="ascii">'
      do n = 1, n_owned
         c = owned_indices(n)
         write(unit_id,'(i0)') c
      end do
      write(unit_id,'(a)') '        </DataArray>'


      ! `fields%mass_flux` is face-centered and owner-to-neighbor oriented.
      ! The volume VTU exports cell-centered rho*u and div(rho*u) diagnostics
      ! so mass-flow information appears in ParaView/PVTU.

      block

         integer :: mf_c, mf_f, mf_owner, mf_nb

         real(rk), allocatable :: mf_div(:)

         allocate(mf_div(mesh%ncells))

         mf_div = 0.0_rk

      

         if (allocated(fields%mass_flux)) then

            do mf_f = 1, mesh%nfaces

               mf_owner = mesh%faces(mf_f)%owner

               mf_nb = mesh%faces(mf_f)%neighbor

               if (mf_nb <= 0) mf_nb = mesh%faces(mf_f)%periodic_neighbor

      

               if (mf_owner > 0) then

                  mf_div(mf_owner) = mf_div(mf_owner) + fields%mass_flux(mf_f) / mesh%cells(mf_owner)%volume

               end if

               if (mf_nb > 0) then

                  mf_div(mf_nb) = mf_div(mf_nb) - fields%mass_flux(mf_f) / mesh%cells(mf_nb)%volume

               end if

            end do

         end if

      

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="mass_flux_vector" NumberOfComponents="3" format="ascii">'

         do mf_c = 1, mesh%ncells

            if (.not. flow%owned(mf_c)) cycle

            write(unit_id,'(3(1x,ES26.16E4))') &

               transport%rho(mf_c) * fields%u(1,mf_c), &

               transport%rho(mf_c) * fields%u(2,mf_c), &

               transport%rho(mf_c) * fields%u(3,mf_c)

         end do

         write(unit_id,'(a)') '        </DataArray>'

      

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="mass_flux_divergence" format="ascii">'

         do mf_c = 1, mesh%ncells

            if (.not. flow%owned(mf_c)) cycle

            write(unit_id,'(1x,ES26.16E4)') mf_div(mf_c)

         end do

         write(unit_id,'(a)') '        </DataArray>'

      

         deallocate(mf_div)

      end block


      ! Current target low-Mach divergence source for div(u)=S projection.

      write(unit_id,'(a)') '        <DataArray type="Float64" Name="lowmach_divergence_source" format="ascii">'

      do c = 1, mesh%ncells

         if (.not. flow%owned(c)) cycle

         if (allocated(fields%divergence_source)) then

            write(unit_id,'(1x,ES26.16E4)') fields%divergence_source(c)

         else

            write(unit_id,'(1x,ES26.16E4)') 0.0_rk

         end if

      end do

      write(unit_id,'(a)') '        </DataArray>'

      if (params%enable_variable_density) then

         call write_lowmach_debug_vtu_arrays(unit_id, mesh, flow, params, fields, transport)

      end if

      write(unit_id,'(a)') '      </CellData>'

      write(unit_id,'(a)') '      <Points>'
      write(unit_id,'(a)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
      do p = 1, mesh%npoints
         write(unit_id,'(3(ES26.16E4,1x))') mesh%points(1,p), mesh%points(2,p), mesh%points(3,p)
      end do
      write(unit_id,'(a)') '        </DataArray>'
      write(unit_id,'(a)') '      </Points>'

      write(unit_id,'(a)') '      <Cells>'

      write(unit_id,'(a)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
      do n = 1, n_owned
         c = owned_indices(n)
         write(unit_id,'(8(i0,1x))') (mesh%cells(c)%nodes(m) - 1, m = 1, 8)
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
      do n = 1, n_owned
         write(unit_id,'(i0)') 8 * n
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a)') '        <DataArray type="UInt8" Name="types" format="ascii">'
      do n = 1, n_owned
         write(unit_id,'(i0)') 12
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a)') '      </Cells>'

      write(unit_id,'(a)') '    </Piece>'
      write(unit_id,'(a)') '  </UnstructuredGrid>'
      write(unit_id,'(a)') '</VTKFile>'

      close(unit_id)

      deallocate(owned_indices)

      ! Rank 0 writes the master PVTU file
            call write_species_energy_conservation_diagnostics(mesh, flow, params, fields, species, energy, transport, &
                                                         step, real(step, rk) * params%dt)
            call write_enthalpy_energy_budget_diagnostics(mesh, flow, params, fields, energy, transport, &
                                                         step, real(step, rk) * params%dt)

if (flow%rank == 0) then
         call write_pvtu_master(params, flow, step)
      end if
   end subroutine write_vtu_unstructured


   !> Writes the master Parallel VTK file (.pvtu) that links the rank pieces.
   subroutine write_pvtu_master(params, flow, step)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow
      integer, intent(in) :: step

      integer :: unit_id, r, k
      character(len=path_len + 64) :: filename
      character(len=64) :: piece_name

      write(filename,'(a,"/VTK/flow_",i6.6,".pvtu")') trim(params%output_dir), step
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a)') '<?xml version="1.0"?>'
      write(unit_id,'(a)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
      write(unit_id,'(a)') '  <PUnstructuredGrid GhostLevel="0">'

      write(unit_id,'(a)') '    <PPointData>'
      write(unit_id,'(a)') '    </PPointData>'

      write(unit_id,'(a)') '    <PCellData Scalars="pressure" Vectors="velocity">'
      write(unit_id,'(a)') '      <PDataArray type="Float64" Name="velocity" NumberOfComponents="3" format="ascii"/>'
      write(unit_id,'(a)') '      <PDataArray type="Float64" Name="pressure" format="ascii"/>'
      write(unit_id,'(a)') '      <PDataArray type="Float64" Name="thermo_pressure" format="ascii"/>'
      write(unit_id,'(a)') '      <PDataArray type="Float64" Name="divergence" format="ascii"/>'
      write(unit_id,'(a)') '      <PDataArray type="Float64" Name="rho" format="ascii"/>'
      write(unit_id,'(a)') '      <PDataArray type="Float64" Name="nu" format="ascii"/>'
      if (params%enable_energy) then
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="temperature" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="enthalpy" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="qrad" format="ascii"/>'
         if (params%enable_species_enthalpy_diffusion) then
            write(unit_id,'(a)') '      <PDataArray type="Float64" Name="species_enthalpy_diffusion" format="ascii"/>'
         end if
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="cp" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="thermal_conductivity" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="thermal_diffusivity" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="rho_thermo" format="ascii"/>'
         if (params%enable_variable_density) then
            write(unit_id,'(a)') '      <PDataArray type="Float64" Name="rho_h_output_state" format="ascii"/>'
            write(unit_id,'(a)') '      <PDataArray type="Float64" Name="rho_h_operator_consistent" format="ascii"/>'
            write(unit_id,'(a)') '      <PDataArray type="Float64" Name="rho_h_density_reconciliation" format="ascii"/>'
            write(unit_id,'(a)') '      <PDataArray type="Float64" Name="relative_rho_h_density_reconciliation" format="ascii"/>'
         end if
      end if
      if (params%enable_species .and. params%nspecies > 0) then
         do k = 1, params%nspecies
            write(unit_id,'(a,a,a)') '      <PDataArray type="Float64" Name="Y_', &
               trim(params%species_name(k)), '" format="ascii"/>'
         end do
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="sum_Y" format="ascii"/>'
         do k = 1, params%nspecies
            write(unit_id,'(a,a,a)') '      <PDataArray type="Float64" Name="D_', &
               trim(params%species_name(k)), '" format="ascii"/>'
         end do
      end if
      write(unit_id,'(a)') '      <PDataArray type="Int32" Name="cell_id" format="ascii"/>'
      write(unit_id,'(a)') '        <PDataArray type="Float64" Name="mass_flux_vector" NumberOfComponents="3"/>'
      write(unit_id,'(a)') '        <PDataArray type="Float64" Name="mass_flux_divergence"/>'
      write(unit_id,'(a)') '        <PDataArray type="Float64" Name="lowmach_divergence_source"/>'

      if (params%enable_variable_density) then
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="lowmach_source_current" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="lowmach_source_projection" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="lowmach_source_difference" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="divu_recomputed" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="divu_minus_S_projection" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="divu_minus_S_current" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="rho_current" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="rho_projection" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="rho_current_minus_projection" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="mass_flux_divergence_recomputed" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="lowmach_source_history_estimate" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="lowmach_source_advective_density" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="u_dot_grad_rho" format="ascii"/>'
         write(unit_id,'(a)') '      <PDataArray type="Float64" Name="continuity_residual_estimate" format="ascii"/>'
      end if
      write(unit_id,'(a)') '    </PCellData>'

      write(unit_id,'(a)') '    <PPoints>'
      write(unit_id,'(a)') '      <PDataArray type="Float64" NumberOfComponents="3" format="ascii"/>'
      write(unit_id,'(a)') '    </PPoints>'

      do r = 0, flow%nprocs - 1
         write(piece_name,'("flow_",i6.6,"_P",i4.4,".vtu")') step, r
         write(unit_id,'(a,a,a)') '    <Piece Source="', trim(piece_name), '"/>'
      end do

      write(unit_id,'(a)') '  </PUnstructuredGrid>'
      write(unit_id,'(a)') '</VTKFile>'

      close(unit_id)
   end subroutine write_pvtu_master


   !> Writes a PVD collection file to allow ParaView to load time-series data.
   subroutine write_pvd_collection(params, flow, nsteps, output_interval, dt)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow
      integer, intent(in) :: nsteps
      integer, intent(in) :: output_interval
      real(rk), intent(in) :: dt

      integer :: unit_id
      integer :: step
      character(len=path_len + 32) :: filename
      character(len=64) :: vtu_name

      if (flow%rank /= 0 .or. .not. params%write_vtu) return

      filename = trim(params%output_dir)//'/VTK/flow.pvd'
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a)') '<?xml version="1.0"?>'
      write(unit_id,'(a)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
      write(unit_id,'(a)') '  <Collection>'

      call write_dataset_line(unit_id, 0, zero)

      do step = output_interval, nsteps, output_interval
         call write_dataset_line(unit_id, step, real(step, rk) * dt)
      end do

      if (mod(nsteps, output_interval) /= 0) then
         call write_dataset_line(unit_id, nsteps, real(nsteps, rk) * dt)
      end if

      write(unit_id,'(a)') '  </Collection>'
      write(unit_id,'(a)') '</VTKFile>'

      close(unit_id)

   contains

      !> Appends a Dataset entry to the PVD collection file.
      !!
      !! @param unit_id The file unit for the .pvd file.
      !! @param step_id The current simulation time step index.
      !! @param time_value The simulation time at this step.
      subroutine write_dataset_line(unit_id, step_id, time_value)
         integer, intent(in) :: unit_id
         integer, intent(in) :: step_id
         real(rk), intent(in) :: time_value

         character(len=32) :: time_text

         write(vtu_name,'("flow_",i6.6,".pvtu")') step_id
         write(time_text,'(ES26.16E4)') time_value
         time_text = adjustl(time_text)

         write(unit_id,'(a,a,a,a,a)') '    <DataSet timestep="', trim(time_text), &
            '" group="" part="0" file="', trim(vtu_name), '"/>'
      end subroutine write_dataset_line

   end subroutine write_pvd_collection


   !> Internal helper to write a scalar field to a VTU file.


   !> Performs sanity checks on hex connectivity before writing output.

   !> Write variable-density low-Mach diagnostics to dedicated CSV files.
   !!
   !! This routine recomputes \( \nabla \cdot u \) from current volumetric face
   !! fluxes and reports source-evolution diagnostics.  The CSV columns named
   !! `*_current` compare against `fields%divergence_source`, which may have
   !! advanced after the projection.  Projection pass/fail should use the
   !! `divu_minus_S_projection_*` diagnostics written by the projection audit
   !! path, where \(S\) is `fields%projection_divergence_source`.
   subroutine write_variable_density_diagnostics(mesh, flow, params, fields, energy, transport, step, time)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: step
      real(rk), intent(in) :: time

      integer :: c, f, lf, fid, owner, nb, ierr, unit_id
      integer :: local_worst_cell, local_winner_rank, global_winner_rank, global_worst_cell
      integer :: mpi_int_send(1), mpi_int_recv(1)
      logical :: file_exists
      logical, save :: vd_diag_initialized = .false.
      logical, save :: vd_worst_initialized = .false.
      character(len=1024) :: filename
      real(rk) :: rho_c, rho_old_c, s_c, div_c, div_res, vol
      real(rk) :: flux_out, mdot_out, mdot_div_c, face_area_value, mass_flux_value
      real(rk) :: local_rho_min, local_rho_max
      real(rk) :: local_rho_old_min, local_rho_old_max
      real(rk) :: local_S_min, local_S_max
      real(rk) :: local_div_min, local_div_max
      real(rk) :: local_res_max, local_res_l2_num
      real(rk) :: local_mdot_div_min, local_mdot_div_max, local_mdot_div_l2_num
      real(rk) :: local_mass, local_net_mdot, local_volume
      real(rk) :: local_h_min, local_h_max
      real(rk) :: local_T_min, local_T_max
      real(rk) :: local_worst_abs, local_worst_divu, local_worst_S
      real(rk) :: local_worst_mdot_div, local_worst_rho, local_worst_T, local_worst_h
      real(rk) :: global_worst_abs, global_worst_divu, global_worst_S
      real(rk) :: global_worst_mdot_div, global_worst_rho, global_worst_T, global_worst_h
      real(rk) :: rho_min, rho_max
      real(rk) :: rho_old_min, rho_old_max
      real(rk) :: S_min, S_max
      real(rk) :: div_min, div_max
      real(rk) :: res_max, res_l2_num, res_l2
      real(rk) :: mdot_div_min, mdot_div_max, mdot_div_l2_num, mdot_div_l2
      real(rk) :: mass_integral, net_boundary_mdot, volume_total
      real(rk) :: h_min, h_max
      real(rk) :: T_min, T_max
      real(rk) :: mpi_reduce_send(1), mpi_reduce_recv(1)

      if (.not. params%enable_variable_density) return
      if (.not. params%write_diagnostics) return
      if (.not. allocated(transport%rho)) return

      local_rho_min = huge(1.0_rk)
      local_rho_max = -huge(1.0_rk)
      local_rho_old_min = huge(1.0_rk)
      local_rho_old_max = -huge(1.0_rk)
      local_S_min = huge(1.0_rk)
      local_S_max = -huge(1.0_rk)
      local_div_min = huge(1.0_rk)
      local_div_max = -huge(1.0_rk)
      local_res_max = 0.0_rk
      local_res_l2_num = 0.0_rk
      local_mdot_div_min = huge(1.0_rk)
      local_mdot_div_max = -huge(1.0_rk)
      local_mdot_div_l2_num = 0.0_rk
      local_mass = 0.0_rk
      local_net_mdot = 0.0_rk
      local_volume = 0.0_rk
      local_h_min = huge(1.0_rk)
      local_h_max = -huge(1.0_rk)
      local_T_min = huge(1.0_rk)
      local_T_max = -huge(1.0_rk)
      local_worst_abs = -1.0_rk
      local_worst_cell = -1
      local_worst_divu = 0.0_rk
      local_worst_S = 0.0_rk
      local_worst_mdot_div = 0.0_rk
      local_worst_rho = 0.0_rk
      local_worst_T = 0.0_rk
      local_worst_h = 0.0_rk

      do c = flow%first_cell, flow%last_cell
         vol = mesh%cells(c)%volume
         rho_c = transport%rho(c)
         rho_old_c = rho_c
         if (allocated(transport%rho_old)) rho_old_c = transport%rho_old(c)

         s_c = 0.0_rk
         if (allocated(fields%divergence_source)) s_c = fields%divergence_source(c)

         div_c = 0.0_rk
         mdot_div_c = 0.0_rk

         if (allocated(fields%face_flux)) then
            do lf = 1, mesh%ncell_faces(c)
               fid = mesh%cell_faces(lf, c)
               if (mesh%faces(fid)%owner == c) then
                  flux_out = fields%face_flux(fid)
               else
                  flux_out = -fields%face_flux(fid)
               end if
               div_c = div_c + flux_out / vol
            end do
         end if

         if (allocated(fields%mass_flux)) then
            do lf = 1, mesh%ncell_faces(c)
               fid = mesh%cell_faces(lf, c)
               if (mesh%faces(fid)%owner == c) then
                  mdot_out = fields%mass_flux(fid)
               else
                  mdot_out = -fields%mass_flux(fid)
               end if
               mdot_div_c = mdot_div_c + mdot_out / vol
            end do
         end if

         div_res = div_c - s_c

         local_rho_min = min(local_rho_min, rho_c)
         local_rho_max = max(local_rho_max, rho_c)
         local_rho_old_min = min(local_rho_old_min, rho_old_c)
         local_rho_old_max = max(local_rho_old_max, rho_old_c)
         local_S_min = min(local_S_min, s_c)
         local_S_max = max(local_S_max, s_c)
         local_div_min = min(local_div_min, div_c)
         local_div_max = max(local_div_max, div_c)
         local_res_max = max(local_res_max, abs(div_res))
         local_res_l2_num = local_res_l2_num + div_res * div_res * vol
         local_mdot_div_min = min(local_mdot_div_min, mdot_div_c)
         local_mdot_div_max = max(local_mdot_div_max, mdot_div_c)
         local_mdot_div_l2_num = local_mdot_div_l2_num + mdot_div_c * mdot_div_c * vol
         local_mass = local_mass + rho_c * vol
         local_volume = local_volume + vol

         if (allocated(energy%h)) then
            local_h_min = min(local_h_min, energy%h(c))
            local_h_max = max(local_h_max, energy%h(c))
         end if
         if (allocated(energy%T)) then
            local_T_min = min(local_T_min, energy%T(c))
            local_T_max = max(local_T_max, energy%T(c))
         end if

         if (abs(div_res) > local_worst_abs) then
            local_worst_abs = abs(div_res)
            local_worst_cell = c
            local_worst_divu = div_c
            local_worst_S = s_c
            local_worst_mdot_div = mdot_div_c
            local_worst_rho = rho_c
            if (allocated(energy%T)) local_worst_T = energy%T(c)
            if (allocated(energy%h)) local_worst_h = energy%h(c)
         end if
      end do

      if (.not. allocated(fields%mass_flux)) then
         local_mdot_div_min = 0.0_rk
         local_mdot_div_max = 0.0_rk
         local_mdot_div_l2_num = 0.0_rk
      end if
      if (.not. allocated(energy%h)) then
         local_h_min = 0.0_rk
         local_h_max = 0.0_rk
      end if
      if (.not. allocated(energy%T)) then
         local_T_min = 0.0_rk
         local_T_max = 0.0_rk
      end if

      if (allocated(fields%mass_flux)) then
         do f = 1, mesh%nfaces
            owner = mesh%faces(f)%owner
            nb = mesh%faces(f)%neighbor
            if (nb <= 0) nb = mesh%faces(f)%periodic_neighbor

            if (nb <= 0) then
               if (owner >= flow%first_cell .and. owner <= flow%last_cell) then
                  local_net_mdot = local_net_mdot + fields%mass_flux(f)
               end if
            end if
         end do
      end if

      mpi_reduce_send(1) = local_rho_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      rho_min = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd rho_min')

      mpi_reduce_send(1) = local_rho_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      rho_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd rho_max')

      mpi_reduce_send(1) = local_rho_old_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      rho_old_min = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd rho_old_min')

      mpi_reduce_send(1) = local_rho_old_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      rho_old_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd rho_old_max')

      mpi_reduce_send(1) = local_S_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      S_min = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd S_min')

      mpi_reduce_send(1) = local_S_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      S_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd S_max')

      mpi_reduce_send(1) = local_div_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      div_min = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd div_min')

      mpi_reduce_send(1) = local_div_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      div_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd div_max')

      mpi_reduce_send(1) = local_res_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      res_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd residual max')

      mpi_reduce_send(1) = local_res_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      res_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd residual l2')

      mpi_reduce_send(1) = local_mdot_div_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      mdot_div_min = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd mdot div min')

      mpi_reduce_send(1) = local_mdot_div_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      mdot_div_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd mdot div max')

      mpi_reduce_send(1) = local_mdot_div_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      mdot_div_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd mdot div l2')

      mpi_reduce_send(1) = local_mass
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      mass_integral = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd mass')

      mpi_reduce_send(1) = local_net_mdot
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      net_boundary_mdot = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd net mdot')

      mpi_reduce_send(1) = local_volume
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      volume_total = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd volume')

      mpi_reduce_send(1) = local_h_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      h_min = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd h_min')

      mpi_reduce_send(1) = local_h_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      h_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd h_max')

      mpi_reduce_send(1) = local_T_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      T_min = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd T_min')

      mpi_reduce_send(1) = local_T_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      T_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd T_max')

      mpi_reduce_send(1) = local_worst_abs
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      global_worst_abs = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd worst residual')

      if (abs(local_worst_abs - global_worst_abs) <= max(1.0e-12_rk, 1.0e-10_rk * max(1.0_rk, global_worst_abs))) then
         local_winner_rank = flow%rank
      else
         local_winner_rank = huge(1)
      end if
      mpi_int_send(1) = local_winner_rank
      call MPI_Allreduce(mpi_int_send, mpi_int_recv, 1, MPI_INTEGER, MPI_MIN, flow%comm, ierr)
      global_winner_rank = mpi_int_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing vd winner rank')

      if (flow%rank == global_winner_rank) then
         mpi_int_send(1) = local_worst_cell
      else
         mpi_int_send(1) = 0
      end if
      call MPI_Allreduce(mpi_int_send, mpi_int_recv, 1, MPI_INTEGER, MPI_SUM, flow%comm, ierr)
      global_worst_cell = mpi_int_recv(1)

      mpi_reduce_send(1) = merge(local_worst_divu, 0.0_rk, flow%rank == global_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_worst_divu = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_worst_S, 0.0_rk, flow%rank == global_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_worst_S = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_worst_mdot_div, 0.0_rk, flow%rank == global_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_worst_mdot_div = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_worst_rho, 0.0_rk, flow%rank == global_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_worst_rho = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_worst_T, 0.0_rk, flow%rank == global_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_worst_T = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_worst_h, 0.0_rk, flow%rank == global_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_worst_h = mpi_reduce_recv(1)

      if (volume_total > 0.0_rk) then
         res_l2 = sqrt(res_l2_num / volume_total)
         mdot_div_l2 = sqrt(mdot_div_l2_num / volume_total)
      else
         res_l2 = 0.0_rk
         mdot_div_l2 = 0.0_rk
      end if

      if (flow%rank == 0) then
         filename = trim(params%output_dir) // '/diagnostics/variable_density_diagnostics.csv'

         if (.not. vd_diag_initialized) then
            file_exists = .false.
            open(newunit=unit_id, file=trim(filename), status='replace', action='write')
            vd_diag_initialized = .true.
         else
            file_exists = .true.
            open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
         end if

         if (.not. file_exists) then
            write(unit_id,'(a)') 'step,time,rho_min,rho_max,rho_old_min,rho_old_max,S_min,S_max,' // &
                                'divu_min,divu_max,divu_minus_S_current_max,divu_minus_S_current_l2,' // &
                                'mass_flux_div_min,mass_flux_div_max,mass_flux_div_l2,' // &
                                'mass_integral,net_boundary_mass_flux,h_min,h_max,T_min,T_max,' // &
                                'worst_residual_abs,worst_residual_rank,worst_residual_cell,' // &
                                'worst_divu,worst_S,worst_mass_flux_div,worst_rho,worst_T,worst_h'
         end if

         write(unit_id,'(i0,21(",",ES26.16E4),2(",",i0),6(",",ES26.16E4))') step, time, rho_min, rho_max, &
              rho_old_min, rho_old_max, S_min, S_max, div_min, div_max, res_max, res_l2, &
              mdot_div_min, mdot_div_max, mdot_div_l2, mass_integral, net_boundary_mdot, &
              h_min, h_max, T_min, T_max, global_worst_abs, global_winner_rank, global_worst_cell, &
              global_worst_divu, global_worst_S, global_worst_mdot_div, global_worst_rho, &
              global_worst_T, global_worst_h
         close(unit_id)
      end if

      if (.not. vd_worst_initialized) then
         if (flow%rank == 0) then
            filename = trim(params%output_dir) // '/diagnostics/variable_density_worst_cell.csv'
            open(newunit=unit_id, file=trim(filename), status='replace', action='write')
            write(unit_id,'(a)') 'step,time,rank,cell,x,y,z,volume,rho,T,h,S,divu,divu_minus_S,mass_flux_div,abs_residual'
            close(unit_id)

            filename = trim(params%output_dir) // '/diagnostics/variable_density_worst_cell_faces.csv'
            open(newunit=unit_id, file=trim(filename), status='replace', action='write')
            write(unit_id,'(a)') 'step,time,rank,cell,local_face,face_id,owner,neighbor,face_area,face_flux,mass_flux,outward_face_flux,outward_mass_flux'
            close(unit_id)
         end if
         vd_worst_initialized = .true.
      end if

      call MPI_Barrier(flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure in vd worst-cell IO barrier')

      if (flow%rank == global_winner_rank .and. local_worst_cell > 0) then
         filename = trim(params%output_dir) // '/diagnostics/variable_density_worst_cell.csv'
         open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
         write(unit_id,'(i0,",",ES26.16E4,",",i0,",",i0,12(",",ES26.16E4))') step, time, flow%rank, local_worst_cell, &
              0.0_rk, 0.0_rk, 0.0_rk, mesh%cells(local_worst_cell)%volume, local_worst_rho, local_worst_T, &
              local_worst_h, local_worst_S, local_worst_divu, local_worst_divu - local_worst_S, &
              local_worst_mdot_div, local_worst_abs
         close(unit_id)

         filename = trim(params%output_dir) // '/diagnostics/variable_density_worst_cell_faces.csv'
         open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
         do lf = 1, mesh%ncell_faces(local_worst_cell)
            fid = mesh%cell_faces(lf, local_worst_cell)
            face_area_value = 0.0_rk

            if (allocated(fields%mass_flux)) then
               mass_flux_value = fields%mass_flux(fid)
            else
               mass_flux_value = 0.0_rk
            end if

            if (mesh%faces(fid)%owner == local_worst_cell) then
               flux_out = fields%face_flux(fid)
               mdot_out = mass_flux_value
            else
               flux_out = -fields%face_flux(fid)
               mdot_out = -mass_flux_value
            end if

            write(unit_id,'(i0,",",ES26.16E4,6(",",i0),5(",",ES26.16E4))') &
                 step, time, flow%rank, local_worst_cell, lf, fid, mesh%faces(fid)%owner, &
                 mesh%faces(fid)%neighbor, face_area_value, fields%face_flux(fid), &
                 mass_flux_value, flux_out, mdot_out
         end do
         close(unit_id)
      end if
      call write_variable_density_boundary_residual_scan(mesh, flow, params, fields, energy, transport, step, time)
      call write_variable_density_projection_audit(mesh, flow, params, fields, energy, transport, step, time)
      call write_variable_density_transport_conservation_diagnostics(mesh, flow, params, fields, energy, transport, step, time)
      call write_variable_density_compatibility_diagnostics(mesh, flow, params, fields, step, time)
      call write_variable_density_continuity_residual_diagnostics(mesh, flow, params, fields, transport, step, time)

   end subroutine write_variable_density_diagnostics

   !> Scan all owned boundary-adjacent cells for low-Mach projection residuals.
   !!
   !! A boundary-adjacent cell is any owned cell touching at least one physical
   !! boundary face, identified by `neighbor <= 0` and `periodic_neighbor <= 0`.
   !! The scan recomputes div(u) and div(rho u) from the current face fluxes and
   !! writes per-rank boundary-cell rows plus a rank-0 global summary.
   subroutine write_variable_density_boundary_residual_scan(mesh, flow, params, fields, energy, transport, step, time)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: step
      real(rk), intent(in) :: time

      integer :: c, lf, fid, nb, ierr, unit_id
      integer :: boundary_face_count
      integer :: local_boundary_count, local_interior_count
      integer :: boundary_count, interior_count
      integer :: local_boundary_worst_cell, local_interior_worst_cell
      integer :: local_boundary_winner_rank, global_boundary_winner_rank
      integer :: local_interior_winner_rank, global_interior_winner_rank
      integer :: global_boundary_worst_cell, global_interior_worst_cell
      integer :: mpi_int_send(1), mpi_int_recv(1)
      logical :: is_boundary_cell
      logical, save :: vd_boundary_scan_initialized = .false.
      character(len=32) :: rank_suffix
      character(len=1024) :: filename
      real(rk) :: vol, rho_c, s_c, div_c, div_res, mdot_div_c
      real(rk) :: flux_out, mdot_out, mass_flux_value
      real(rk) :: T_c, h_c
      real(rk) :: local_boundary_volume, local_interior_volume
      real(rk) :: boundary_volume, interior_volume
      real(rk) :: local_boundary_res_max, local_boundary_res_l2_num
      real(rk) :: local_interior_res_max, local_interior_res_l2_num
      real(rk) :: boundary_res_max, boundary_res_l2_num, boundary_res_l2
      real(rk) :: interior_res_max, interior_res_l2_num, interior_res_l2
      real(rk) :: local_boundary_div_min, local_boundary_div_max
      real(rk) :: local_boundary_S_min, local_boundary_S_max
      real(rk) :: local_boundary_mdot_min, local_boundary_mdot_max
      real(rk) :: boundary_div_min, boundary_div_max
      real(rk) :: boundary_S_min, boundary_S_max
      real(rk) :: boundary_mdot_min, boundary_mdot_max
      real(rk) :: local_boundary_worst_divu, local_boundary_worst_S
      real(rk) :: local_boundary_worst_mdot_div, local_boundary_worst_rho
      real(rk) :: local_boundary_worst_T, local_boundary_worst_h
      real(rk) :: local_interior_worst_divu, local_interior_worst_S
      real(rk) :: local_interior_worst_mdot_div, local_interior_worst_rho
      real(rk) :: local_interior_worst_T, local_interior_worst_h
      real(rk) :: global_boundary_worst_divu, global_boundary_worst_S
      real(rk) :: global_boundary_worst_mdot_div, global_boundary_worst_rho
      real(rk) :: global_boundary_worst_T, global_boundary_worst_h
      real(rk) :: global_interior_worst_divu, global_interior_worst_S
      real(rk) :: global_interior_worst_mdot_div, global_interior_worst_rho
      real(rk) :: global_interior_worst_T, global_interior_worst_h
      real(rk) :: mpi_reduce_send(1), mpi_reduce_recv(1)

      if (.not. params%enable_variable_density) return
      if (.not. params%write_diagnostics) return
      if (.not. allocated(transport%rho)) return

      write(rank_suffix, '(i0)') flow%rank

      filename = trim(params%output_dir) // '/diagnostics/variable_density_boundary_residual_cells_rank' // &
                 trim(adjustl(rank_suffix)) // '.csv'
      if (.not. vd_boundary_scan_initialized) then
         open(newunit=unit_id, file=trim(filename), status='replace', action='write')
         write(unit_id,'(a)') 'step,time,rank,cell,boundary_face_count,volume,rho,T,h,S,divu,divu_minus_S,mass_flux_div,abs_residual'
      else
         open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
      end if

      local_boundary_count = 0
      local_interior_count = 0
      local_boundary_volume = 0.0_rk
      local_interior_volume = 0.0_rk
      local_boundary_res_max = 0.0_rk
      local_interior_res_max = 0.0_rk
      local_boundary_res_l2_num = 0.0_rk
      local_interior_res_l2_num = 0.0_rk

      local_boundary_div_min = huge(1.0_rk)
      local_boundary_div_max = -huge(1.0_rk)
      local_boundary_S_min = huge(1.0_rk)
      local_boundary_S_max = -huge(1.0_rk)
      local_boundary_mdot_min = huge(1.0_rk)
      local_boundary_mdot_max = -huge(1.0_rk)

      local_boundary_worst_cell = -1
      local_interior_worst_cell = -1
      local_boundary_worst_divu = 0.0_rk
      local_boundary_worst_S = 0.0_rk
      local_boundary_worst_mdot_div = 0.0_rk
      local_boundary_worst_rho = 0.0_rk
      local_boundary_worst_T = 0.0_rk
      local_boundary_worst_h = 0.0_rk
      local_interior_worst_divu = 0.0_rk
      local_interior_worst_S = 0.0_rk
      local_interior_worst_mdot_div = 0.0_rk
      local_interior_worst_rho = 0.0_rk
      local_interior_worst_T = 0.0_rk
      local_interior_worst_h = 0.0_rk

      do c = flow%first_cell, flow%last_cell
         vol = mesh%cells(c)%volume
         rho_c = transport%rho(c)

         s_c = 0.0_rk
         if (allocated(fields%divergence_source)) s_c = fields%divergence_source(c)

         T_c = 0.0_rk
         h_c = 0.0_rk
         if (allocated(energy%T)) T_c = energy%T(c)
         if (allocated(energy%h)) h_c = energy%h(c)

         div_c = 0.0_rk
         mdot_div_c = 0.0_rk
         boundary_face_count = 0

         do lf = 1, mesh%ncell_faces(c)
            fid = mesh%cell_faces(lf, c)

            nb = mesh%faces(fid)%neighbor
            if (nb <= 0) then
               if (mesh%faces(fid)%periodic_neighbor <= 0) boundary_face_count = boundary_face_count + 1
            end if

            if (allocated(fields%face_flux)) then
               if (mesh%faces(fid)%owner == c) then
                  flux_out = fields%face_flux(fid)
               else
                  flux_out = -fields%face_flux(fid)
               end if
               div_c = div_c + flux_out / vol
            end if

            if (allocated(fields%mass_flux)) then
               mass_flux_value = fields%mass_flux(fid)
               if (mesh%faces(fid)%owner == c) then
                  mdot_out = mass_flux_value
               else
                  mdot_out = -mass_flux_value
               end if
               mdot_div_c = mdot_div_c + mdot_out / vol
            end if
         end do

         is_boundary_cell = boundary_face_count > 0
         div_res = div_c - s_c

         if (is_boundary_cell) then
            local_boundary_count = local_boundary_count + 1
            local_boundary_volume = local_boundary_volume + vol
            local_boundary_res_l2_num = local_boundary_res_l2_num + div_res * div_res * vol
            local_boundary_div_min = min(local_boundary_div_min, div_c)
            local_boundary_div_max = max(local_boundary_div_max, div_c)
            local_boundary_S_min = min(local_boundary_S_min, s_c)
            local_boundary_S_max = max(local_boundary_S_max, s_c)
            local_boundary_mdot_min = min(local_boundary_mdot_min, mdot_div_c)
            local_boundary_mdot_max = max(local_boundary_mdot_max, mdot_div_c)

            write(unit_id,'(i0,",",ES26.16E4,",",i0,",",i0,",",i0,9(",",ES26.16E4))') &
                 step, time, flow%rank, c, boundary_face_count, vol, rho_c, T_c, h_c, s_c, &
                 div_c, div_res, mdot_div_c, abs(div_res)

            if (abs(div_res) > local_boundary_res_max) then
               local_boundary_res_max = abs(div_res)
               local_boundary_worst_cell = c
               local_boundary_worst_divu = div_c
               local_boundary_worst_S = s_c
               local_boundary_worst_mdot_div = mdot_div_c
               local_boundary_worst_rho = rho_c
               local_boundary_worst_T = T_c
               local_boundary_worst_h = h_c
            end if
         else
            local_interior_count = local_interior_count + 1
            local_interior_volume = local_interior_volume + vol
            local_interior_res_l2_num = local_interior_res_l2_num + div_res * div_res * vol

            if (abs(div_res) > local_interior_res_max) then
               local_interior_res_max = abs(div_res)
               local_interior_worst_cell = c
               local_interior_worst_divu = div_c
               local_interior_worst_S = s_c
               local_interior_worst_mdot_div = mdot_div_c
               local_interior_worst_rho = rho_c
               local_interior_worst_T = T_c
               local_interior_worst_h = h_c
            end if
         end if
      end do

      close(unit_id)

      if (local_boundary_count == 0) then
         local_boundary_div_min = 0.0_rk
         local_boundary_div_max = 0.0_rk
         local_boundary_S_min = 0.0_rk
         local_boundary_S_max = 0.0_rk
         local_boundary_mdot_min = 0.0_rk
         local_boundary_mdot_max = 0.0_rk
      end if

      mpi_int_send(1) = local_boundary_count
      call MPI_Allreduce(mpi_int_send, mpi_int_recv, 1, MPI_INTEGER, MPI_SUM, flow%comm, ierr)
      boundary_count = mpi_int_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing boundary cell count')

      mpi_int_send(1) = local_interior_count
      call MPI_Allreduce(mpi_int_send, mpi_int_recv, 1, MPI_INTEGER, MPI_SUM, flow%comm, ierr)
      interior_count = mpi_int_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing interior cell count')

      mpi_reduce_send(1) = local_boundary_volume
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      boundary_volume = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_interior_volume
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      interior_volume = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_boundary_res_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      boundary_res_max = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_interior_res_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      interior_res_max = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_boundary_res_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      boundary_res_l2_num = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_interior_res_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      interior_res_l2_num = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_boundary_div_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      boundary_div_min = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_boundary_div_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      boundary_div_max = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_boundary_S_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      boundary_S_min = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_boundary_S_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      boundary_S_max = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_boundary_mdot_min
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      boundary_mdot_min = mpi_reduce_recv(1)

      mpi_reduce_send(1) = local_boundary_mdot_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      boundary_mdot_max = mpi_reduce_recv(1)

      if (boundary_volume > 0.0_rk) then
         boundary_res_l2 = sqrt(boundary_res_l2_num / boundary_volume)
      else
         boundary_res_l2 = 0.0_rk
      end if

      if (interior_volume > 0.0_rk) then
         interior_res_l2 = sqrt(interior_res_l2_num / interior_volume)
      else
         interior_res_l2 = 0.0_rk
      end if

      if (abs(local_boundary_res_max - boundary_res_max) <= max(1.0e-12_rk, 1.0e-10_rk * max(1.0_rk, boundary_res_max))) then
         local_boundary_winner_rank = flow%rank
      else
         local_boundary_winner_rank = huge(1)
      end if

      mpi_int_send(1) = local_boundary_winner_rank
      call MPI_Allreduce(mpi_int_send, mpi_int_recv, 1, MPI_INTEGER, MPI_MIN, flow%comm, ierr)
      global_boundary_winner_rank = mpi_int_recv(1)

      if (abs(local_interior_res_max - interior_res_max) <= max(1.0e-12_rk, 1.0e-10_rk * max(1.0_rk, interior_res_max))) then
         local_interior_winner_rank = flow%rank
      else
         local_interior_winner_rank = huge(1)
      end if

      mpi_int_send(1) = local_interior_winner_rank
      call MPI_Allreduce(mpi_int_send, mpi_int_recv, 1, MPI_INTEGER, MPI_MIN, flow%comm, ierr)
      global_interior_winner_rank = mpi_int_recv(1)

      if (flow%rank == global_boundary_winner_rank) then
         mpi_int_send(1) = local_boundary_worst_cell
      else
         mpi_int_send(1) = 0
      end if
      call MPI_Allreduce(mpi_int_send, mpi_int_recv, 1, MPI_INTEGER, MPI_SUM, flow%comm, ierr)
      global_boundary_worst_cell = mpi_int_recv(1)

      if (flow%rank == global_interior_winner_rank) then
         mpi_int_send(1) = local_interior_worst_cell
      else
         mpi_int_send(1) = 0
      end if
      call MPI_Allreduce(mpi_int_send, mpi_int_recv, 1, MPI_INTEGER, MPI_SUM, flow%comm, ierr)
      global_interior_worst_cell = mpi_int_recv(1)

      mpi_reduce_send(1) = merge(local_boundary_worst_divu, 0.0_rk, flow%rank == global_boundary_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_boundary_worst_divu = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_boundary_worst_S, 0.0_rk, flow%rank == global_boundary_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_boundary_worst_S = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_boundary_worst_mdot_div, 0.0_rk, flow%rank == global_boundary_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_boundary_worst_mdot_div = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_boundary_worst_rho, 0.0_rk, flow%rank == global_boundary_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_boundary_worst_rho = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_boundary_worst_T, 0.0_rk, flow%rank == global_boundary_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_boundary_worst_T = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_boundary_worst_h, 0.0_rk, flow%rank == global_boundary_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_boundary_worst_h = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_interior_worst_divu, 0.0_rk, flow%rank == global_interior_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_interior_worst_divu = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_interior_worst_S, 0.0_rk, flow%rank == global_interior_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_interior_worst_S = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_interior_worst_mdot_div, 0.0_rk, flow%rank == global_interior_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_interior_worst_mdot_div = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_interior_worst_rho, 0.0_rk, flow%rank == global_interior_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_interior_worst_rho = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_interior_worst_T, 0.0_rk, flow%rank == global_interior_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_interior_worst_T = mpi_reduce_recv(1)

      mpi_reduce_send(1) = merge(local_interior_worst_h, 0.0_rk, flow%rank == global_interior_winner_rank)
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      global_interior_worst_h = mpi_reduce_recv(1)

      if (flow%rank == 0) then
         filename = trim(params%output_dir) // '/diagnostics/variable_density_boundary_residual_summary.csv'
         if (.not. vd_boundary_scan_initialized) then
            open(newunit=unit_id, file=trim(filename), status='replace', action='write')
            write(unit_id,'(a)') 'step,time,boundary_cell_count,interior_cell_count,boundary_res_max,boundary_res_l2,interior_res_max,interior_res_l2,boundary_divu_min,boundary_divu_max,boundary_S_min,boundary_S_max,boundary_mass_flux_div_min,boundary_mass_flux_div_max,boundary_worst_rank,boundary_worst_cell,boundary_worst_divu,boundary_worst_S,boundary_worst_mass_flux_div,boundary_worst_rho,boundary_worst_T,boundary_worst_h,interior_worst_rank,interior_worst_cell,interior_worst_divu,interior_worst_S,interior_worst_mass_flux_div,interior_worst_rho,interior_worst_T,interior_worst_h'
         else
            open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
         end if

         write(unit_id,'(i0,",",ES26.16E4,2(",",i0),10(",",ES26.16E4),2(",",i0),6(",",ES26.16E4),2(",",i0),6(",",ES26.16E4))') &
              step, time, boundary_count, interior_count, boundary_res_max, boundary_res_l2, &
              interior_res_max, interior_res_l2, boundary_div_min, boundary_div_max, boundary_S_min, boundary_S_max, &
              boundary_mdot_min, boundary_mdot_max, global_boundary_winner_rank, global_boundary_worst_cell, &
              global_boundary_worst_divu, global_boundary_worst_S, global_boundary_worst_mdot_div, &
              global_boundary_worst_rho, global_boundary_worst_T, global_boundary_worst_h, &
              global_interior_winner_rank, global_interior_worst_cell, global_interior_worst_divu, &
              global_interior_worst_S, global_interior_worst_mdot_div, global_interior_worst_rho, &
              global_interior_worst_T, global_interior_worst_h
         close(unit_id)
      end if

      vd_boundary_scan_initialized = .true.

   end subroutine write_variable_density_boundary_residual_scan

   !> Targeted projection audit for variable-density low-Mach residuals.
   !!
   !! Writes the local worst residual cell on each rank plus cell 1 when owned.
   !! This is intended to diagnose corner/boundary pressure-correction problems
   !! without changing the numerical scheme.
   subroutine write_variable_density_projection_audit(mesh, flow, params, fields, energy, transport, step, time)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: step
      real(rk), intent(in) :: time

      integer :: c, unit_cells, unit_faces
      integer :: local_worst_cell, boundary_face_count, physical_boundary_face_count
      logical, save :: projection_audit_initialized = .false.
      character(len=32) :: rank_suffix
      character(len=1024) :: cells_file, faces_file
      real(rk) :: local_worst_abs
      real(rk) :: s_c, div_c, div_res, mdot_div_c

      if (.not. params%enable_variable_density) return
      if (.not. params%write_diagnostics) return
      if (.not. allocated(transport%rho)) return

      local_worst_abs = -1.0_rk
      local_worst_cell = -1

      do c = flow%first_cell, flow%last_cell
         call compute_projection_audit_cell_balance(mesh, fields, c, div_c, mdot_div_c, &
                                                    boundary_face_count, physical_boundary_face_count)

         s_c = 0.0_rk
         if (allocated(fields%divergence_source)) s_c = fields%divergence_source(c)

         div_res = div_c - s_c
         if (abs(div_res) > local_worst_abs) then
            local_worst_abs = abs(div_res)
            local_worst_cell = c
         end if
      end do

      write(rank_suffix, '(i0)') flow%rank
      cells_file = trim(params%output_dir) // '/diagnostics/variable_density_projection_audit_cells_rank' // &
                   trim(adjustl(rank_suffix)) // '.csv'
      faces_file = trim(params%output_dir) // '/diagnostics/variable_density_projection_audit_faces_rank' // &
                   trim(adjustl(rank_suffix)) // '.csv'

      if (.not. projection_audit_initialized) then
         open(newunit=unit_cells, file=trim(cells_file), status='replace', action='write')
         write(unit_cells,'(a)') 'step,time,rank,audit_label,cell,volume,rho,T,h,S,divu,divu_minus_S,' // &
                                 'mass_flux_div,required_flux_sum,actual_flux_sum,flux_defect,' // &
                                 'mass_flux_sum,boundary_face_count,physical_boundary_face_count,pressure'
         close(unit_cells)

         open(newunit=unit_faces, file=trim(faces_file), status='replace', action='write')
         write(unit_faces,'(a)') 'step,time,rank,audit_label,cell,local_face,face_id,owner,neighbor,periodic_neighbor,' // &
                                 'is_boundary_face,is_physical_boundary_face,owner_pressure,cell_pressure,neighbor_pressure,' // &
                                 'dp_neighbor_minus_cell,face_flux,outward_face_flux,mass_flux,outward_mass_flux'
         close(unit_faces)
         projection_audit_initialized = .true.
      end if

      open(newunit=unit_cells, file=trim(cells_file), status='unknown', position='append', action='write')
      open(newunit=unit_faces, file=trim(faces_file), status='unknown', position='append', action='write')

      if (local_worst_cell > 0) then
         call write_projection_audit_one_cell(mesh, flow, fields, energy, transport, step, time, &
                                              local_worst_cell, 'local_worst', unit_cells, unit_faces)
      end if

      if (flow%first_cell <= 1 .and. 1 <= flow%last_cell) then
         if (local_worst_cell /= 1) then
            call write_projection_audit_one_cell(mesh, flow, fields, energy, transport, step, time, &
                                                 1, 'cell_1', unit_cells, unit_faces)
         end if
      end if

      close(unit_cells)
      close(unit_faces)

   end subroutine write_variable_density_projection_audit


   !> Compute recomputed divergence and mass-flux divergence for one cell.
   subroutine compute_projection_audit_cell_balance(mesh, fields, c, div_c, mdot_div_c, &
                                                    boundary_face_count, physical_boundary_face_count)
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(in) :: fields
      integer, intent(in) :: c
      real(rk), intent(out) :: div_c
      real(rk), intent(out) :: mdot_div_c
      integer, intent(out) :: boundary_face_count
      integer, intent(out) :: physical_boundary_face_count

      integer :: lf, fid, nb
      real(rk) :: vol, flux_out, mdot_out, mass_flux_value

      div_c = 0.0_rk
      mdot_div_c = 0.0_rk
      boundary_face_count = 0
      physical_boundary_face_count = 0

      vol = mesh%cells(c)%volume
      if (vol <= 0.0_rk) return

      do lf = 1, mesh%ncell_faces(c)
         fid = mesh%cell_faces(lf, c)

         nb = mesh%faces(fid)%neighbor
         if (nb <= 0) then
            boundary_face_count = boundary_face_count + 1
            if (mesh%faces(fid)%periodic_neighbor <= 0) physical_boundary_face_count = physical_boundary_face_count + 1
         end if

         if (allocated(fields%face_flux)) then
            if (mesh%faces(fid)%owner == c) then
               flux_out = fields%face_flux(fid)
            else
               flux_out = -fields%face_flux(fid)
            end if
            div_c = div_c + flux_out / vol
         end if

         if (allocated(fields%mass_flux)) then
            mass_flux_value = fields%mass_flux(fid)
            if (mesh%faces(fid)%owner == c) then
               mdot_out = mass_flux_value
            else
               mdot_out = -mass_flux_value
            end if
            mdot_div_c = mdot_div_c + mdot_out / vol
         end if
      end do

   end subroutine compute_projection_audit_cell_balance


   !> Write one audited cell and all its faces.
   subroutine write_projection_audit_one_cell(mesh, flow, fields, energy, transport, step, time, &
                                              c, audit_label, unit_cells, unit_faces)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(flow_fields_t), intent(in) :: fields
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: step
      real(rk), intent(in) :: time
      integer, intent(in) :: c
      character(len=*), intent(in) :: audit_label
      integer, intent(in) :: unit_cells
      integer, intent(in) :: unit_faces

      integer :: lf, fid, nb, periodic_nb
      integer :: boundary_face_count, physical_boundary_face_count
      integer :: is_boundary_face, is_physical_boundary_face
      real(rk) :: vol, rho_c, T_c, h_c, s_c, div_c, div_res, mdot_div_c
      real(rk) :: required_flux_sum, actual_flux_sum, flux_defect, mass_flux_sum
      real(rk) :: cell_pressure, owner_pressure, neighbor_pressure, dp_neighbor_minus_cell
      real(rk) :: face_flux_value, outward_face_flux, mass_flux_value, outward_mass_flux

      if (c < 1) return

      vol = mesh%cells(c)%volume
      rho_c = transport%rho(c)
      T_c = 0.0_rk
      h_c = 0.0_rk
      s_c = 0.0_rk

      if (allocated(energy%T)) T_c = energy%T(c)
      if (allocated(energy%h)) h_c = energy%h(c)
      if (allocated(fields%divergence_source)) s_c = fields%divergence_source(c)

      call compute_projection_audit_cell_balance(mesh, fields, c, div_c, mdot_div_c, &
                                                 boundary_face_count, physical_boundary_face_count)

      div_res = div_c - s_c
      required_flux_sum = s_c * vol
      actual_flux_sum = div_c * vol
      flux_defect = actual_flux_sum - required_flux_sum
      mass_flux_sum = mdot_div_c * vol
      cell_pressure = merge(fields%p(c), 0.0_rk, allocated(fields%p))

      write(unit_cells,'(i0,",",ES26.16E4,",",i0,",",a,",",i0,12(",",ES26.16E4),2(",",i0),",",ES26.16E4)') &
           step, time, flow%rank, trim(audit_label), c, vol, rho_c, T_c, h_c, s_c, div_c, div_res, &
           mdot_div_c, required_flux_sum, actual_flux_sum, flux_defect, mass_flux_sum, &
           boundary_face_count, physical_boundary_face_count, cell_pressure

      do lf = 1, mesh%ncell_faces(c)
         fid = mesh%cell_faces(lf, c)
         nb = mesh%faces(fid)%neighbor
         periodic_nb = mesh%faces(fid)%periodic_neighbor

         is_boundary_face = 0
         is_physical_boundary_face = 0
         if (nb <= 0) then
            is_boundary_face = 1
            if (periodic_nb <= 0) is_physical_boundary_face = 1
         end if

         face_flux_value = 0.0_rk
         if (allocated(fields%face_flux)) face_flux_value = fields%face_flux(fid)

         mass_flux_value = 0.0_rk
         if (allocated(fields%mass_flux)) mass_flux_value = fields%mass_flux(fid)

         if (mesh%faces(fid)%owner == c) then
            outward_face_flux = face_flux_value
            outward_mass_flux = mass_flux_value
         else
            outward_face_flux = -face_flux_value
            outward_mass_flux = -mass_flux_value
         end if

         owner_pressure = merge(fields%p(mesh%faces(fid)%owner), 0.0_rk, allocated(fields%p) .and. mesh%faces(fid)%owner > 0)
         neighbor_pressure = merge(fields%p(nb), 0.0_rk, allocated(fields%p) .and. nb > 0)
         dp_neighbor_minus_cell = neighbor_pressure - cell_pressure

         write(unit_faces,'(i0,",",ES26.16E4,",",i0,",",a,8(",",i0),8(",",ES26.16E4))') &
              step, time, flow%rank, trim(audit_label), c, lf, fid, mesh%faces(fid)%owner, nb, periodic_nb, &
              is_boundary_face, is_physical_boundary_face, owner_pressure, cell_pressure, neighbor_pressure, &
              dp_neighbor_minus_cell, face_flux_value, outward_face_flux, mass_flux_value, outward_mass_flux
      end do

   end subroutine write_projection_audit_one_cell

   !> Write low-Mach compatibility and source-time-level diagnostics.
   !!
   !! The current divergence_source may be advanced by the energy/thermo density
   !! sync after the projection has already used the previous source level.  The
   !! projection_divergence_source field stores the source that was actually
   !! consumed by the latest projection RHS and outlet balancing.  Comparing both
   !! residuals distinguishes pressure/projection error from source time-level
   !! mismatch.
   !! CSV columns with `_current` compare against `fields%divergence_source`
   !! after the energy/thermo update; columns with `_projection` compare
   !! against the source copied before projection.
   subroutine write_variable_density_compatibility_diagnostics(mesh, flow, params, fields, step, time)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      integer, intent(in) :: step
      real(rk), intent(in) :: time

      integer :: c, f, lf, fid, owner, nb, ierr, unit_id
      logical, save :: compatibility_diag_initialized = .false.
      character(len=1024) :: filename
      real(rk) :: vol, s_c, s_proj_c, div_c, res_current, res_projection, s_delta, flux_out
      real(rk) :: local_integral_S_dV, integral_S_dV
      real(rk) :: local_integral_S_projection_dV, integral_S_projection_dV
      real(rk) :: local_net_boundary_volume_flux, net_boundary_volume_flux
      real(rk) :: compatibility_error, compatibility_error_projection
      real(rk) :: local_volume, volume_total, mean_S, mean_S_projection
      real(rk) :: local_S_l2_num, S_l2_num, S_l2
      real(rk) :: local_S_projection_l2_num, S_projection_l2_num, S_projection_l2
      real(rk) :: local_S_abs_max, S_abs_max
      real(rk) :: local_S_projection_abs_max, S_projection_abs_max
      real(rk) :: local_res_max, res_max
      real(rk) :: local_res_projection_max, res_projection_max
      real(rk) :: local_res_l2_num, res_l2_num, res_l2
      real(rk) :: local_res_projection_l2_num, res_projection_l2_num, res_projection_l2
      real(rk) :: local_S_delta_max, S_delta_max
      real(rk) :: local_S_delta_l2_num, S_delta_l2_num, S_delta_l2
      real(rk) :: rel_res_max, rel_res_l2
      real(rk) :: rel_res_projection_max, rel_res_projection_l2
      real(rk) :: denom
      real(rk) :: mpi_reduce_send(1), mpi_reduce_recv(1)

      if (.not. params%enable_variable_density) return
      if (.not. params%write_diagnostics) return

      local_integral_S_dV = 0.0_rk
      local_integral_S_projection_dV = 0.0_rk
      local_net_boundary_volume_flux = 0.0_rk
      local_volume = 0.0_rk
      local_S_l2_num = 0.0_rk
      local_S_projection_l2_num = 0.0_rk
      local_S_abs_max = 0.0_rk
      local_S_projection_abs_max = 0.0_rk
      local_res_max = 0.0_rk
      local_res_projection_max = 0.0_rk
      local_res_l2_num = 0.0_rk
      local_res_projection_l2_num = 0.0_rk
      local_S_delta_max = 0.0_rk
      local_S_delta_l2_num = 0.0_rk

      do c = flow%first_cell, flow%last_cell
         vol = mesh%cells(c)%volume

         s_c = 0.0_rk
         if (allocated(fields%divergence_source)) s_c = fields%divergence_source(c)

         s_proj_c = s_c
         if (allocated(fields%projection_divergence_source)) s_proj_c = fields%projection_divergence_source(c)

         div_c = 0.0_rk
         if (allocated(fields%face_flux)) then
            do lf = 1, mesh%ncell_faces(c)
               fid = mesh%cell_faces(lf, c)
               if (mesh%faces(fid)%owner == c) then
                  flux_out = fields%face_flux(fid)
               else
                  flux_out = -fields%face_flux(fid)
               end if
               div_c = div_c + flux_out / vol
            end do
         end if

         res_current = div_c - s_c
         res_projection = div_c - s_proj_c
         s_delta = s_c - s_proj_c

         local_integral_S_dV = local_integral_S_dV + s_c * vol
         local_integral_S_projection_dV = local_integral_S_projection_dV + s_proj_c * vol
         local_volume = local_volume + vol

         local_S_l2_num = local_S_l2_num + s_c * s_c * vol
         local_S_projection_l2_num = local_S_projection_l2_num + s_proj_c * s_proj_c * vol
         local_S_abs_max = max(local_S_abs_max, abs(s_c))
         local_S_projection_abs_max = max(local_S_projection_abs_max, abs(s_proj_c))

         local_res_max = max(local_res_max, abs(res_current))
         local_res_projection_max = max(local_res_projection_max, abs(res_projection))
         local_res_l2_num = local_res_l2_num + res_current * res_current * vol
         local_res_projection_l2_num = local_res_projection_l2_num + res_projection * res_projection * vol

         local_S_delta_max = max(local_S_delta_max, abs(s_delta))
         local_S_delta_l2_num = local_S_delta_l2_num + s_delta * s_delta * vol
      end do

      if (allocated(fields%face_flux)) then
         do f = 1, mesh%nfaces
            owner = mesh%faces(f)%owner
            nb = mesh%faces(f)%neighbor
            if (nb <= 0) nb = mesh%faces(f)%periodic_neighbor

            if (nb <= 0) then
               if (owner >= flow%first_cell .and. owner <= flow%last_cell) then
                  local_net_boundary_volume_flux = local_net_boundary_volume_flux + fields%face_flux(f)
               end if
            end if
         end do
      end if

      mpi_reduce_send(1) = local_integral_S_dV
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      integral_S_dV = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility integral S dV')

      mpi_reduce_send(1) = local_integral_S_projection_dV
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      integral_S_projection_dV = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility projection integral S dV')

      mpi_reduce_send(1) = local_net_boundary_volume_flux
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      net_boundary_volume_flux = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility boundary flux')

      mpi_reduce_send(1) = local_volume
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      volume_total = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility volume')

      mpi_reduce_send(1) = local_S_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      S_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility S l2')

      mpi_reduce_send(1) = local_S_projection_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      S_projection_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility projection S l2')

      mpi_reduce_send(1) = local_S_abs_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      S_abs_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility S max')

      mpi_reduce_send(1) = local_S_projection_abs_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      S_projection_abs_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility projection S max')

      mpi_reduce_send(1) = local_res_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      res_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility residual max')

      mpi_reduce_send(1) = local_res_projection_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      res_projection_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility projection residual max')

      mpi_reduce_send(1) = local_res_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      res_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility residual l2')

      mpi_reduce_send(1) = local_res_projection_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      res_projection_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility projection residual l2')

      mpi_reduce_send(1) = local_S_delta_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      S_delta_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility source delta max')

      mpi_reduce_send(1) = local_S_delta_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      S_delta_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing compatibility source delta l2')

      compatibility_error = net_boundary_volume_flux - integral_S_dV
      compatibility_error_projection = net_boundary_volume_flux - integral_S_projection_dV

      if (volume_total > 0.0_rk) then
         mean_S = integral_S_dV / volume_total
         mean_S_projection = integral_S_projection_dV / volume_total
         S_l2 = sqrt(S_l2_num / volume_total)
         S_projection_l2 = sqrt(S_projection_l2_num / volume_total)
         res_l2 = sqrt(res_l2_num / volume_total)
         res_projection_l2 = sqrt(res_projection_l2_num / volume_total)
         S_delta_l2 = sqrt(S_delta_l2_num / volume_total)
      else
         mean_S = 0.0_rk
         mean_S_projection = 0.0_rk
         S_l2 = 0.0_rk
         S_projection_l2 = 0.0_rk
         res_l2 = 0.0_rk
         res_projection_l2 = 0.0_rk
         S_delta_l2 = 0.0_rk
      end if

      denom = max(S_abs_max, tiny(1.0_rk))
      rel_res_max = res_max / denom

      denom = max(S_l2, tiny(1.0_rk))
      rel_res_l2 = res_l2 / denom

      denom = max(S_projection_abs_max, tiny(1.0_rk))
      rel_res_projection_max = res_projection_max / denom

      denom = max(S_projection_l2, tiny(1.0_rk))
      rel_res_projection_l2 = res_projection_l2 / denom

      if (flow%rank == 0) then
         filename = trim(params%output_dir) // '/diagnostics/variable_density_compatibility.csv'

         if (.not. compatibility_diag_initialized) then
            open(newunit=unit_id, file=trim(filename), status='replace', action='write')
            compatibility_diag_initialized = .true.
            write(unit_id,'(a)') 'step,time,integral_S_current_dV,net_boundary_volume_flux,' // &
                                'net_boundary_volume_flux_minus_integral_S_current_dV,volume_total,mean_S_current,' // &
                                'S_current_l2,S_current_abs_max,divu_minus_S_current_max,divu_minus_S_current_l2,' // &
                                'relative_divu_minus_S_current_max,relative_divu_minus_S_current_l2,' // &
                                'integral_S_projection_dV,' // &
                                'net_boundary_volume_flux_minus_integral_S_projection_dV,' // &
                                'mean_S_projection,S_projection_l2,S_projection_abs_max,' // &
                                'divu_minus_S_projection_max,divu_minus_S_projection_l2,' // &
                                'relative_divu_minus_S_projection_max,' // &
                                'relative_divu_minus_S_projection_l2,' // &
                                'S_current_minus_S_projection_max,S_current_minus_S_projection_l2'
         else
            open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
         end if

         write(unit_id,'(i0,23(",",ES26.16E4))') step, time, integral_S_dV, net_boundary_volume_flux, &
                                              compatibility_error, volume_total, mean_S, S_l2, S_abs_max, &
                                              res_max, res_l2, rel_res_max, rel_res_l2, &
                                              integral_S_projection_dV, compatibility_error_projection, &
                                              mean_S_projection, S_projection_l2, S_projection_abs_max, &
                                              res_projection_max, res_projection_l2, &
                                              rel_res_projection_max, rel_res_projection_l2, &
                                              S_delta_max, S_delta_l2
         close(unit_id)
      end if

   end subroutine write_variable_density_compatibility_diagnostics

   !> Write variable-density transport/conservation diagnostics.
   !!
   !! Boundary mass flux is positive outward, so global conservative closure is:
   !!
   !!   d/dt integral(rho dV) + net_boundary_mass_flux = 0
   !!
   !! The time derivative here is output-to-output, not per-step.  This keeps
   !! the diagnostic independent of the integrator state.
   !> Write variable-density mass/transport conservation diagnostics with
   !! explicit density and mass-flux time levels.
   !!
   !! Sign convention:
   !!   positive boundary flux is outward from the owner cell/domain.
   !!
   !! For a conservative domain mass balance:
   !!   dM/dt + net_boundary_mass_flux = 0
   !!
   !! This routine reports several versions of the boundary mass flux so the
   !! diagnostics can distinguish a true conservative transport error from a
   !! current-density/projection-density time-level mismatch.
   subroutine write_variable_density_transport_conservation_diagnostics(mesh, flow, params, fields, energy, transport, step, time)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: step
      real(rk), intent(in) :: time

      integer :: c, f, owner, nb, ierr, unit_id
      logical, save :: transport_diag_initialized = .false.
      character(len=1024) :: filename
      real(rk), save :: prev_time = 0.0_rk
      real(rk), save :: prev_mass_current = 0.0_rk
      real(rk), save :: prev_mass_projection = 0.0_rk
      real(rk), save :: prev_rho_h_current = 0.0_rk
      real(rk) :: vol, rho_c, rho_p, h_c, T_c
      real(rk) :: face_flux, mass_flux_stored
      real(rk) :: dt_prev
      real(rk) :: delta_mass_current, mass_rate_current
      real(rk) :: delta_mass_projection, mass_rate_projection
      real(rk) :: delta_rho_h_current, rho_h_rate_current
      real(rk) :: defect_current_stored, defect_current_current_rho
      real(rk) :: defect_projection_projection_rho
      real(rk) :: rel_defect_current_stored, rel_defect_current_current_rho
      real(rk) :: rel_defect_projection_projection_rho
      real(rk) :: denom
      real(rk) :: sum_send(12), sum_recv(12)
      real(rk) :: min_send(4), min_recv(4)
      real(rk) :: max_send(5), max_recv(5)
      real(rk) :: local_volume
      real(rk) :: local_mass_current, local_mass_projection
      real(rk) :: local_net_boundary_volume_flux
      real(rk) :: local_net_boundary_mass_flux_stored
      real(rk) :: local_net_boundary_mass_flux_current_rho
      real(rk) :: local_net_boundary_mass_flux_projection_rho
      real(rk) :: local_rho_current_l2_num, local_rho_projection_l2_num
      real(rk) :: local_rho_diff_l2_num
      real(rk) :: local_rho_h_current, local_rho_h_projection
      real(rk) :: volume_total
      real(rk) :: mass_current, mass_projection
      real(rk) :: current_minus_projection_mass
      real(rk) :: net_boundary_volume_flux
      real(rk) :: net_boundary_mass_flux_stored
      real(rk) :: net_boundary_mass_flux_current_rho
      real(rk) :: net_boundary_mass_flux_projection_rho
      real(rk) :: stored_minus_current_rho_boundary_mass_flux
      real(rk) :: stored_minus_projection_rho_boundary_mass_flux
      real(rk) :: rho_current_min, rho_current_max, rho_current_mean, rho_current_l2
      real(rk) :: rho_projection_min, rho_projection_max, rho_projection_mean, rho_projection_l2
      real(rk) :: rho_diff_max, rho_diff_l2
      real(rk) :: h_min, h_max, T_min, T_max
      real(rk) :: rho_h_current, rho_h_projection

      if (.not. params%enable_variable_density) return
      if (.not. params%write_diagnostics) return
      if (.not. allocated(transport%rho)) return

      local_volume = 0.0_rk
      local_mass_current = 0.0_rk
      local_mass_projection = 0.0_rk
      local_net_boundary_volume_flux = 0.0_rk
      local_net_boundary_mass_flux_stored = 0.0_rk
      local_net_boundary_mass_flux_current_rho = 0.0_rk
      local_net_boundary_mass_flux_projection_rho = 0.0_rk
      local_rho_current_l2_num = 0.0_rk
      local_rho_projection_l2_num = 0.0_rk
      local_rho_diff_l2_num = 0.0_rk
      local_rho_h_current = 0.0_rk
      local_rho_h_projection = 0.0_rk

      rho_current_min = huge(1.0_rk)
      rho_projection_min = huge(1.0_rk)
      h_min = huge(1.0_rk)
      T_min = huge(1.0_rk)
      rho_current_max = -huge(1.0_rk)
      rho_projection_max = -huge(1.0_rk)
      rho_diff_max = 0.0_rk
      h_max = -huge(1.0_rk)
      T_max = -huge(1.0_rk)

      do c = flow%first_cell, flow%last_cell
         vol = mesh%cells(c)%volume
         rho_c = transport%rho(c)
         rho_p = rho_c
         if (allocated(fields%projection_rho)) then
            if (fields%projection_rho(c) > 0.0_rk) rho_p = fields%projection_rho(c)
         end if

         h_c = 0.0_rk
         if (allocated(energy%h)) h_c = energy%h(c)
         T_c = 0.0_rk
         if (allocated(energy%T)) T_c = energy%T(c)

         local_volume = local_volume + vol
         local_mass_current = local_mass_current + rho_c * vol
         local_mass_projection = local_mass_projection + rho_p * vol
         local_rho_current_l2_num = local_rho_current_l2_num + rho_c * rho_c * vol
         local_rho_projection_l2_num = local_rho_projection_l2_num + rho_p * rho_p * vol
         local_rho_diff_l2_num = local_rho_diff_l2_num + (rho_c - rho_p) * (rho_c - rho_p) * vol
         local_rho_h_current = local_rho_h_current + rho_c * h_c * vol
         local_rho_h_projection = local_rho_h_projection + rho_p * h_c * vol

         rho_current_min = min(rho_current_min, rho_c)
         rho_current_max = max(rho_current_max, rho_c)
         rho_projection_min = min(rho_projection_min, rho_p)
         rho_projection_max = max(rho_projection_max, rho_p)
         rho_diff_max = max(rho_diff_max, abs(rho_c - rho_p))
         h_min = min(h_min, h_c)
         h_max = max(h_max, h_c)
         T_min = min(T_min, T_c)
         T_max = max(T_max, T_c)
      end do

      if (allocated(fields%face_flux)) then
         do f = 1, mesh%nfaces
            owner = mesh%faces(f)%owner
            nb = mesh%faces(f)%neighbor
            if (nb <= 0) nb = mesh%faces(f)%periodic_neighbor

            if (nb <= 0) then
               if (owner >= flow%first_cell .and. owner <= flow%last_cell) then
                  face_flux = fields%face_flux(f)
                  rho_c = transport%rho(owner)
                  rho_p = rho_c
                  if (allocated(fields%projection_rho)) then
                     if (fields%projection_rho(owner) > 0.0_rk) rho_p = fields%projection_rho(owner)
                  end if

                  local_net_boundary_volume_flux = local_net_boundary_volume_flux + face_flux
                  local_net_boundary_mass_flux_current_rho = local_net_boundary_mass_flux_current_rho + face_flux * rho_c
                  local_net_boundary_mass_flux_projection_rho = local_net_boundary_mass_flux_projection_rho + face_flux * rho_p

                  if (allocated(fields%mass_flux)) then
                     mass_flux_stored = fields%mass_flux(f)
                     local_net_boundary_mass_flux_stored = local_net_boundary_mass_flux_stored + mass_flux_stored
                  end if
               end if
            end if
         end do
      end if

      sum_send = 0.0_rk
      sum_send(1) = local_volume
      sum_send(2) = local_mass_current
      sum_send(3) = local_mass_projection
      sum_send(4) = local_net_boundary_volume_flux
      sum_send(5) = local_net_boundary_mass_flux_stored
      sum_send(6) = local_net_boundary_mass_flux_current_rho
      sum_send(7) = local_net_boundary_mass_flux_projection_rho
      sum_send(8) = local_rho_current_l2_num
      sum_send(9) = local_rho_projection_l2_num
      sum_send(10) = local_rho_diff_l2_num
      sum_send(11) = local_rho_h_current
      sum_send(12) = local_rho_h_projection
      call MPI_Allreduce(sum_send, sum_recv, 12, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing transport conservation sums')

      min_send(1) = rho_current_min
      min_send(2) = rho_projection_min
      min_send(3) = h_min
      min_send(4) = T_min
      call MPI_Allreduce(min_send, min_recv, 4, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing transport conservation minima')

      max_send(1) = rho_current_max
      max_send(2) = rho_projection_max
      max_send(3) = rho_diff_max
      max_send(4) = h_max
      max_send(5) = T_max
      call MPI_Allreduce(max_send, max_recv, 5, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing transport conservation maxima')

      volume_total = sum_recv(1)
      mass_current = sum_recv(2)
      mass_projection = sum_recv(3)
      current_minus_projection_mass = mass_current - mass_projection
      net_boundary_volume_flux = sum_recv(4)
      net_boundary_mass_flux_stored = sum_recv(5)
      net_boundary_mass_flux_current_rho = sum_recv(6)
      net_boundary_mass_flux_projection_rho = sum_recv(7)
      stored_minus_current_rho_boundary_mass_flux = net_boundary_mass_flux_stored - net_boundary_mass_flux_current_rho
      stored_minus_projection_rho_boundary_mass_flux = net_boundary_mass_flux_stored - net_boundary_mass_flux_projection_rho
      rho_h_current = sum_recv(11)
      rho_h_projection = sum_recv(12)

      if (volume_total > 0.0_rk) then
         rho_current_mean = mass_current / volume_total
         rho_projection_mean = mass_projection / volume_total
         rho_current_l2 = sqrt(sum_recv(8) / volume_total)
         rho_projection_l2 = sqrt(sum_recv(9) / volume_total)
         rho_diff_l2 = sqrt(sum_recv(10) / volume_total)
      else
         rho_current_mean = 0.0_rk
         rho_projection_mean = 0.0_rk
         rho_current_l2 = 0.0_rk
         rho_projection_l2 = 0.0_rk
         rho_diff_l2 = 0.0_rk
      end if

      rho_current_min = min_recv(1)
      rho_projection_min = min_recv(2)
      h_min = min_recv(3)
      T_min = min_recv(4)
      rho_current_max = max_recv(1)
      rho_projection_max = max_recv(2)
      rho_diff_max = max_recv(3)
      h_max = max_recv(4)
      T_max = max_recv(5)

      if (transport_diag_initialized) then
         dt_prev = time - prev_time
      else
         dt_prev = 0.0_rk
      end if

      if (dt_prev > tiny(1.0_rk)) then
         delta_mass_current = mass_current - prev_mass_current
         mass_rate_current = delta_mass_current / dt_prev
         delta_mass_projection = mass_projection - prev_mass_projection
         mass_rate_projection = delta_mass_projection / dt_prev
         delta_rho_h_current = rho_h_current - prev_rho_h_current
         rho_h_rate_current = delta_rho_h_current / dt_prev
      else
         delta_mass_current = 0.0_rk
         mass_rate_current = 0.0_rk
         delta_mass_projection = 0.0_rk
         mass_rate_projection = 0.0_rk
         delta_rho_h_current = 0.0_rk
         rho_h_rate_current = 0.0_rk
      end if

      defect_current_stored = mass_rate_current + net_boundary_mass_flux_stored
      defect_current_current_rho = mass_rate_current + net_boundary_mass_flux_current_rho
      defect_projection_projection_rho = mass_rate_projection + net_boundary_mass_flux_projection_rho

      denom = max(max(abs(mass_rate_current), abs(net_boundary_mass_flux_stored)), tiny(1.0_rk))
      rel_defect_current_stored = abs(defect_current_stored) / denom
      denom = max(max(abs(mass_rate_current), abs(net_boundary_mass_flux_current_rho)), tiny(1.0_rk))
      rel_defect_current_current_rho = abs(defect_current_current_rho) / denom
      denom = max(max(abs(mass_rate_projection), abs(net_boundary_mass_flux_projection_rho)), tiny(1.0_rk))
      rel_defect_projection_projection_rho = abs(defect_projection_projection_rho) / denom

      if (flow%rank == 0) then
         filename = trim(params%output_dir) // '/diagnostics/variable_density_transport_conservation.csv'

         if (.not. transport_diag_initialized) then
            open(newunit=unit_id, file=trim(filename), status='replace', action='write')
            write(unit_id,'(a)') 'step,time,volume_total,' // &
                                'mass_integral_current,mass_integral_projection,current_minus_projection_mass_integral,' // &
                                'net_boundary_volume_flux,net_boundary_mass_flux_stored,' // &
                                'net_boundary_mass_flux_current_rho,net_boundary_mass_flux_projection_rho,' // &
                                'stored_minus_current_rho_boundary_mass_flux,stored_minus_projection_rho_boundary_mass_flux,' // &
                                'delta_time_since_previous,delta_mass_current_since_previous,mass_rate_current_since_previous,' // &
                                'mass_balance_defect_current_stored_flux,relative_mass_balance_defect_current_stored_flux,' // &
                                'mass_balance_defect_current_current_rho_flux,relative_mass_balance_defect_current_current_rho_flux,' // &
                                'delta_mass_projection_since_previous,mass_rate_projection_since_previous,' // &
                                'mass_balance_defect_projection_projection_rho_flux,' // &
                                'relative_mass_balance_defect_projection_projection_rho_flux,' // &
                                'rho_current_min,rho_current_max,rho_current_mean,rho_current_l2,' // &
                                'rho_projection_min,rho_projection_max,rho_projection_mean,rho_projection_l2,' // &
                                'rho_current_minus_projection_max,rho_current_minus_projection_l2,' // &
                                'h_min,h_max,T_min,T_max,' // &
                                'rho_h_integral_current,rho_h_integral_projection,' // &
                                'delta_rho_h_current_since_previous,rho_h_current_rate_since_previous'
         else
            open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
         end if

         write(unit_id,'(i0,40(",",ES26.16E4))') step, time, volume_total, &
            mass_current, mass_projection, current_minus_projection_mass, &
            net_boundary_volume_flux, net_boundary_mass_flux_stored, &
            net_boundary_mass_flux_current_rho, net_boundary_mass_flux_projection_rho, &
            stored_minus_current_rho_boundary_mass_flux, stored_minus_projection_rho_boundary_mass_flux, &
            dt_prev, delta_mass_current, mass_rate_current, &
            defect_current_stored, rel_defect_current_stored, &
            defect_current_current_rho, rel_defect_current_current_rho, &
            delta_mass_projection, mass_rate_projection, &
            defect_projection_projection_rho, rel_defect_projection_projection_rho, &
            rho_current_min, rho_current_max, rho_current_mean, rho_current_l2, &
            rho_projection_min, rho_projection_max, rho_projection_mean, rho_projection_l2, &
            rho_diff_max, rho_diff_l2, h_min, h_max, T_min, T_max, &
            rho_h_current, rho_h_projection, delta_rho_h_current, rho_h_rate_current
         close(unit_id)
      end if

      transport_diag_initialized = .true.
      prev_time = time
      prev_mass_current = mass_current
      prev_mass_projection = mass_projection
      prev_rho_h_current = rho_h_current

   end subroutine write_variable_density_transport_conservation_diagnostics

   !> Write conservative continuity residual diagnostics for variable-density mode.
   !!
   !! This diagnostics-only routine separates the local density-history
   !! source relation from conservative mass continuity:
   !!
   !!   d(rho)/dt + div(rho u) = 0
   !!
   !! and its expanded finite-volume form:
   !!
   !!   d(rho)/dt + rho div(u) + u.grad(rho) = 0
   !!
   !! with u.grad(rho) estimated as div(rho u) - rho div(u).
   subroutine write_variable_density_continuity_residual_diagnostics(mesh, flow, params, fields, transport, step, time)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: step
      real(rk), intent(in) :: time

      integer :: c, lf, fid, ierr, unit_id
      logical, save :: continuity_diag_initialized = .false.
      character(len=1024) :: filename
      real(rk) :: vol, rho_current, rho_projection
      real(rk) :: flux_out, mass_flux_out
      real(rk) :: divu, div_mass_flux, drho_dt
      real(rk) :: rho_divu, udotgradrho
      real(rk) :: history_res, conservative_res, expanded_res
      real(rk) :: local_volume, volume_total
      real(rk) :: local_integral_drho_dt_dV, integral_drho_dt_dV
      real(rk) :: local_integral_rho_divu_dV, integral_rho_divu_dV
      real(rk) :: local_integral_udotgradrho_dV, integral_udotgradrho_dV
      real(rk) :: local_integral_div_mass_flux_dV, integral_div_mass_flux_dV
      real(rk) :: local_integral_history_res_dV, integral_history_res_dV
      real(rk) :: local_integral_conservative_res_dV, integral_conservative_res_dV
      real(rk) :: local_integral_expanded_res_dV, integral_expanded_res_dV
      real(rk) :: local_history_res_max, history_res_max
      real(rk) :: local_conservative_res_max, conservative_res_max
      real(rk) :: local_expanded_res_max, expanded_res_max
      real(rk) :: local_udotgradrho_max, udotgradrho_max
      real(rk) :: local_history_res_l2_num, history_res_l2_num, history_res_l2
      real(rk) :: local_conservative_res_l2_num, conservative_res_l2_num, conservative_res_l2
      real(rk) :: local_expanded_res_l2_num, expanded_res_l2_num, expanded_res_l2
      real(rk) :: local_udotgradrho_l2_num, udotgradrho_l2_num, udotgradrho_l2
      real(rk) :: local_drho_dt_l2_num, drho_dt_l2_num, drho_dt_l2
      real(rk) :: local_div_mass_flux_l2_num, div_mass_flux_l2_num, div_mass_flux_l2
      real(rk) :: local_rho_divu_l2_num, rho_divu_l2_num, rho_divu_l2
      real(rk) :: denom, rel_conservative_res_l2, rel_conservative_res_max
      real(rk) :: mpi_reduce_send(1), mpi_reduce_recv(1)

      if (.not. params%enable_variable_density) return
      if (.not. params%write_diagnostics) return
      if (.not. allocated(transport%rho)) return
      if (.not. allocated(fields%face_flux)) return

      local_volume = 0.0_rk
      local_integral_drho_dt_dV = 0.0_rk
      local_integral_rho_divu_dV = 0.0_rk
      local_integral_udotgradrho_dV = 0.0_rk
      local_integral_div_mass_flux_dV = 0.0_rk
      local_integral_history_res_dV = 0.0_rk
      local_integral_conservative_res_dV = 0.0_rk
      local_integral_expanded_res_dV = 0.0_rk
      local_history_res_max = 0.0_rk
      local_conservative_res_max = 0.0_rk
      local_expanded_res_max = 0.0_rk
      local_udotgradrho_max = 0.0_rk
      local_history_res_l2_num = 0.0_rk
      local_conservative_res_l2_num = 0.0_rk
      local_expanded_res_l2_num = 0.0_rk
      local_udotgradrho_l2_num = 0.0_rk
      local_drho_dt_l2_num = 0.0_rk
      local_div_mass_flux_l2_num = 0.0_rk
      local_rho_divu_l2_num = 0.0_rk

      do c = flow%first_cell, flow%last_cell
         vol = mesh%cells(c)%volume
         rho_current = transport%rho(c)
         rho_projection = rho_current
         if (allocated(fields%projection_rho)) rho_projection = fields%projection_rho(c)

         divu = 0.0_rk
         div_mass_flux = 0.0_rk

         do lf = 1, mesh%ncell_faces(c)
            fid = mesh%cell_faces(lf, c)
            if (mesh%faces(fid)%owner == c) then
               flux_out = fields%face_flux(fid)
               if (allocated(fields%mass_flux)) then
                  mass_flux_out = fields%mass_flux(fid)
               else
                  mass_flux_out = rho_projection * fields%face_flux(fid)
               end if
            else
               flux_out = -fields%face_flux(fid)
               if (allocated(fields%mass_flux)) then
                  mass_flux_out = -fields%mass_flux(fid)
               else
                  mass_flux_out = -rho_projection * fields%face_flux(fid)
               end if
            end if
            divu = divu + flux_out / vol
            div_mass_flux = div_mass_flux + mass_flux_out / vol
         end do

         if (params%dt > 0.0_rk) then
            drho_dt = (rho_current - rho_projection) / params%dt
         else
            drho_dt = 0.0_rk
         end if

         rho_divu = rho_current * divu
         udotgradrho = div_mass_flux - rho_divu

         history_res = drho_dt + rho_divu
         conservative_res = drho_dt + div_mass_flux
         expanded_res = drho_dt + rho_divu + udotgradrho

         local_volume = local_volume + vol
         local_integral_drho_dt_dV = local_integral_drho_dt_dV + drho_dt * vol
         local_integral_rho_divu_dV = local_integral_rho_divu_dV + rho_divu * vol
         local_integral_udotgradrho_dV = local_integral_udotgradrho_dV + udotgradrho * vol
         local_integral_div_mass_flux_dV = local_integral_div_mass_flux_dV + div_mass_flux * vol
         local_integral_history_res_dV = local_integral_history_res_dV + history_res * vol
         local_integral_conservative_res_dV = local_integral_conservative_res_dV + conservative_res * vol
         local_integral_expanded_res_dV = local_integral_expanded_res_dV + expanded_res * vol

         local_history_res_max = max(local_history_res_max, abs(history_res))
         local_conservative_res_max = max(local_conservative_res_max, abs(conservative_res))
         local_expanded_res_max = max(local_expanded_res_max, abs(expanded_res))
         local_udotgradrho_max = max(local_udotgradrho_max, abs(udotgradrho))

         local_history_res_l2_num = local_history_res_l2_num + history_res * history_res * vol
         local_conservative_res_l2_num = local_conservative_res_l2_num + conservative_res * conservative_res * vol
         local_expanded_res_l2_num = local_expanded_res_l2_num + expanded_res * expanded_res * vol
         local_udotgradrho_l2_num = local_udotgradrho_l2_num + udotgradrho * udotgradrho * vol
         local_drho_dt_l2_num = local_drho_dt_l2_num + drho_dt * drho_dt * vol
         local_div_mass_flux_l2_num = local_div_mass_flux_l2_num + div_mass_flux * div_mass_flux * vol
         local_rho_divu_l2_num = local_rho_divu_l2_num + rho_divu * rho_divu * vol
      end do

      mpi_reduce_send(1) = local_volume
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      volume_total = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity volume')

      mpi_reduce_send(1) = local_integral_drho_dt_dV
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      integral_drho_dt_dV = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity drho_dt')

      mpi_reduce_send(1) = local_integral_rho_divu_dV
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      integral_rho_divu_dV = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity rho_divu')

      mpi_reduce_send(1) = local_integral_udotgradrho_dV
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      integral_udotgradrho_dV = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity udotgradrho')

      mpi_reduce_send(1) = local_integral_div_mass_flux_dV
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      integral_div_mass_flux_dV = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity div_mass_flux')

      mpi_reduce_send(1) = local_integral_history_res_dV
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      integral_history_res_dV = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity history residual')

      mpi_reduce_send(1) = local_integral_conservative_res_dV
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      integral_conservative_res_dV = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity conservative residual')

      mpi_reduce_send(1) = local_integral_expanded_res_dV
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      integral_expanded_res_dV = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity expanded residual')

      mpi_reduce_send(1) = local_history_res_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      history_res_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity history max')

      mpi_reduce_send(1) = local_conservative_res_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      conservative_res_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity conservative max')

      mpi_reduce_send(1) = local_expanded_res_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      expanded_res_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity expanded max')

      mpi_reduce_send(1) = local_udotgradrho_max
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      udotgradrho_max = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity udotgradrho max')

      mpi_reduce_send(1) = local_history_res_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      history_res_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity history l2')

      mpi_reduce_send(1) = local_conservative_res_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      conservative_res_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity conservative l2')

      mpi_reduce_send(1) = local_expanded_res_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      expanded_res_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity expanded l2')

      mpi_reduce_send(1) = local_udotgradrho_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      udotgradrho_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity udotgradrho l2')

      mpi_reduce_send(1) = local_drho_dt_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      drho_dt_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity drho_dt l2')

      mpi_reduce_send(1) = local_div_mass_flux_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      div_mass_flux_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity div_mass_flux l2')

      mpi_reduce_send(1) = local_rho_divu_l2_num
      call MPI_Allreduce(mpi_reduce_send, mpi_reduce_recv, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      rho_divu_l2_num = mpi_reduce_recv(1)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing continuity rho_divu l2')

      if (volume_total > 0.0_rk) then
         history_res_l2 = sqrt(history_res_l2_num / volume_total)
         conservative_res_l2 = sqrt(conservative_res_l2_num / volume_total)
         expanded_res_l2 = sqrt(expanded_res_l2_num / volume_total)
         udotgradrho_l2 = sqrt(udotgradrho_l2_num / volume_total)
         drho_dt_l2 = sqrt(drho_dt_l2_num / volume_total)
         div_mass_flux_l2 = sqrt(div_mass_flux_l2_num / volume_total)
         rho_divu_l2 = sqrt(rho_divu_l2_num / volume_total)
      else
         history_res_l2 = 0.0_rk
         conservative_res_l2 = 0.0_rk
         expanded_res_l2 = 0.0_rk
         udotgradrho_l2 = 0.0_rk
         drho_dt_l2 = 0.0_rk
         div_mass_flux_l2 = 0.0_rk
         rho_divu_l2 = 0.0_rk
      end if

      denom = max(div_mass_flux_l2 + drho_dt_l2, tiny(1.0_rk))
      rel_conservative_res_l2 = conservative_res_l2 / denom

      denom = max(max(abs(integral_div_mass_flux_dV), abs(integral_drho_dt_dV)), tiny(1.0_rk))
      rel_conservative_res_max = abs(integral_conservative_res_dV) / denom

      if (flow%rank == 0) then
         filename = trim(params%output_dir) // '/diagnostics/variable_density_continuity_residual.csv'

         if (.not. continuity_diag_initialized) then
            open(newunit=unit_id, file=trim(filename), status='replace', action='write')
            continuity_diag_initialized = .true.
            write(unit_id,'(a)') 'step,time,volume_total,' // &
                                'integral_drho_dt_dV,integral_rho_current_divu_dV,' // &
                                'integral_udotgradrho_dV,integral_div_mass_flux_dV,' // &
                                'integral_drho_dt_plus_rho_current_divu_dV,' // &
                                'integral_drho_dt_plus_div_mass_flux_dV,' // &
                                'integral_drho_dt_plus_rho_current_divu_plus_udotgradrho_dV,' // &
                                'history_residual_max,history_residual_l2,' // &
                                'conservative_residual_max,conservative_residual_l2,' // &
                                'expanded_residual_max,expanded_residual_l2,' // &
                                'udotgradrho_max,udotgradrho_l2,' // &
                                'drho_dt_l2,div_mass_flux_l2,rho_current_divu_l2,' // &
                                'relative_conservative_residual_l2,' // &
                                'relative_integral_conservative_residual'
         else
            open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
         end if

         write(unit_id,'(i0,22(",",ES26.16E4))') step, time, volume_total, &
            integral_drho_dt_dV, integral_rho_divu_dV, integral_udotgradrho_dV, &
            integral_div_mass_flux_dV, integral_history_res_dV, &
            integral_conservative_res_dV, integral_expanded_res_dV, &
            history_res_max, history_res_l2, conservative_res_max, conservative_res_l2, &
            expanded_res_max, expanded_res_l2, udotgradrho_max, udotgradrho_l2, &
            drho_dt_l2, div_mass_flux_l2, rho_divu_l2, rel_conservative_res_l2, &
            rel_conservative_res_max
         close(unit_id)
      end if

   end subroutine write_variable_density_continuity_residual_diagnostics


   !> Append local rho*h density-reconciliation fields to VTU CellData.
   !!
   !! These fields are spatial diagnostics for the global energy-density
   !! reconciliation.  They are not global closure metrics; those remain in
   !! diagnostics/enthalpy_energy_budget.csv.
   subroutine write_energy_reconciliation_vtu_arrays(unit_id, mesh, flow, energy, transport)
      integer, intent(in) :: unit_id
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport

      call write_energy_reconciliation_scalar(unit_id, mesh, flow, energy, transport, &
                                              'rho_h_output_state')
      call write_energy_reconciliation_scalar(unit_id, mesh, flow, energy, transport, &
                                              'rho_h_operator_consistent')
      call write_energy_reconciliation_scalar(unit_id, mesh, flow, energy, transport, &
                                              'rho_h_density_reconciliation')
      call write_energy_reconciliation_scalar(unit_id, mesh, flow, energy, transport, &
                                              'relative_rho_h_density_reconciliation')

   end subroutine write_energy_reconciliation_vtu_arrays


   subroutine write_energy_reconciliation_scalar(unit_id, mesh, flow, energy, transport, name)
      integer, intent(in) :: unit_id
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      character(len=*), intent(in) :: name

      integer :: c
      real(rk) :: value

      write(unit_id,'(a)') '        <DataArray type="Float64" Name="' // trim(name) // '" format="ascii">'
      do c = 1, mesh%ncells
         if (.not. flow%owned(c)) cycle
         value = energy_reconciliation_value(energy, transport, c, trim(name))
         write(unit_id,'(ES26.16E4)') value
      end do
      write(unit_id,'(a)') '        </DataArray>'

   end subroutine write_energy_reconciliation_scalar


   real(rk) function energy_reconciliation_value(energy, transport, c, name) result(value)
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: c
      character(len=*), intent(in) :: name

      real(rk) :: rho_h_output, rho_h_operator, rho_h_recon, denom

      value = 0.0_rk
      rho_h_output = 0.0_rk
      rho_h_operator = 0.0_rk
      rho_h_recon = 0.0_rk

      if (.not. allocated(energy%h)) return
      if (.not. allocated(transport%rho)) return
      if (c < 1 .or. c > size(energy%h) .or. c > size(transport%rho)) return

      rho_h_output = transport%rho(c) * energy%h(c)

      ! Before the first energy update, no operator-consistent cellwise state
      ! exists.  Fall back to the output-state value so step-0 visualization has
      ! a defined, zero reconciliation field.
      rho_h_operator = rho_h_output
      if (energy%operator_consistent_rho_h_available == 1 .and. &
          allocated(energy%operator_consistent_rho_h)) then
         if (size(energy%operator_consistent_rho_h) >= c) then
            rho_h_operator = energy%operator_consistent_rho_h(c)
         end if
      end if

      rho_h_recon = rho_h_output - rho_h_operator

      select case (trim(name))
      case ('rho_h_output_state')
         value = rho_h_output
      case ('rho_h_operator_consistent')
         value = rho_h_operator
      case ('rho_h_density_reconciliation')
         value = rho_h_recon
      case ('relative_rho_h_density_reconciliation')
         denom = max(abs(rho_h_output), abs(rho_h_operator), tiny(1.0_rk))
         value = abs(rho_h_recon) / denom
      case default
         value = 0.0_rk
      end select

   end function energy_reconciliation_value

   !> Append variable-density low-Mach debug fields to the VTU CellData block.
   !!
   !! This helper writes one scalar value per owned cell, matching the normal
   !! parallel VTU partitioning.  It intentionally does not change numerics.
   subroutine write_lowmach_debug_vtu_arrays(unit_id, mesh, flow, params, fields, transport)
      integer, intent(in) :: unit_id
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(transport_properties_t), intent(in) :: transport

      if (.not. params%enable_variable_density) return

      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'lowmach_source_current')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'lowmach_source_projection')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'lowmach_source_difference')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'divu_recomputed')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'divu_minus_S_projection')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'divu_minus_S_current')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'rho_current')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'rho_projection')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'rho_current_minus_projection')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'mass_flux_divergence_recomputed')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'lowmach_source_history_estimate')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'lowmach_source_advective_density')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'u_dot_grad_rho')
      call write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, 'continuity_residual_estimate')

   end subroutine write_lowmach_debug_vtu_arrays


   subroutine write_lowmach_debug_scalar(unit_id, mesh, flow, fields, transport, name)
      integer, intent(in) :: unit_id
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(flow_fields_t), intent(in) :: fields
      type(transport_properties_t), intent(in) :: transport
      character(len=*), intent(in) :: name

      integer :: c
      real(rk) :: value

      write(unit_id,'(a)') '        <DataArray type="Float64" Name="' // trim(name) // '" format="ascii">'
      do c = 1, mesh%ncells
         if (.not. flow%owned(c)) cycle
         value = lowmach_debug_value(mesh, fields, transport, c, trim(name))
         write(unit_id,'(es24.16)') value
      end do
      write(unit_id,'(a)') '        </DataArray>'

   end subroutine write_lowmach_debug_scalar


   real(rk) function lowmach_debug_value(mesh, fields, transport, c, name) result(value)
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(in) :: fields
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: c
      character(len=*), intent(in) :: name

      real(rk) :: divu, s_current, s_projection
      real(rk) :: rho_current, rho_projection, rho_floor
      real(rk) :: div_mass, ugrad_rho, source_advective, source_history
      real(rk) :: continuity_residual

      value = 0.0_rk
      rho_floor = tiny(1.0_rk)

      divu = lowmach_cell_divu(mesh, fields, c)
      div_mass = lowmach_cell_div_mass_flux(mesh, fields, c)
      s_current = lowmach_cell_source_current(fields, c)
      s_projection = lowmach_cell_source_projection(fields, c)
      rho_current = lowmach_cell_rho_current(transport, c)
      rho_projection = lowmach_cell_rho_projection(fields, transport, c)

      ugrad_rho = lowmach_cell_u_dot_grad_rho(mesh, fields, transport, c)
      if (rho_current > rho_floor) then
         source_advective = -ugrad_rho / rho_current
      else
         source_advective = 0.0_rk
      end if
      source_history = s_current - source_advective

      ! Since rho_old is advanced after source construction, this visualization
      ! estimate reconstructs d(rho)/dt from the history-source component:
      ! S_history = (rho_old - rho)/(rho dt) = -drho_dt/rho.
      continuity_residual = -rho_current * source_history + div_mass

      select case (trim(name))
      case ('lowmach_source_current')
         value = s_current
      case ('lowmach_source_projection')
         value = s_projection
      case ('lowmach_source_difference')
         value = s_current - s_projection
      case ('divu_recomputed')
         value = divu
      case ('divu_minus_S_projection')
         value = divu - s_projection
      case ('divu_minus_S_current')
         value = divu - s_current
      case ('rho_current')
         value = rho_current
      case ('rho_projection')
         value = rho_projection
      case ('rho_current_minus_projection')
         value = rho_current - rho_projection
      case ('mass_flux_divergence_recomputed')
         value = div_mass
      case ('lowmach_source_history_estimate')
         value = source_history
      case ('lowmach_source_advective_density')
         value = source_advective
      case ('u_dot_grad_rho')
         value = ugrad_rho
      case ('continuity_residual_estimate')
         value = continuity_residual
      case default
         value = 0.0_rk
      end select

   end function lowmach_debug_value


   real(rk) function lowmach_cell_source_current(fields, c) result(value)
      type(flow_fields_t), intent(in) :: fields
      integer, intent(in) :: c

      value = 0.0_rk
      if (allocated(fields%divergence_source)) then
         if (size(fields%divergence_source) >= c) value = fields%divergence_source(c)
      end if
   end function lowmach_cell_source_current


   real(rk) function lowmach_cell_source_projection(fields, c) result(value)
      type(flow_fields_t), intent(in) :: fields
      integer, intent(in) :: c

      value = lowmach_cell_source_current(fields, c)
      if (allocated(fields%projection_divergence_source)) then
         if (size(fields%projection_divergence_source) >= c) value = fields%projection_divergence_source(c)
      end if
   end function lowmach_cell_source_projection


   real(rk) function lowmach_cell_rho_current(transport, c) result(value)
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: c

      value = 0.0_rk
      if (allocated(transport%rho)) then
         if (size(transport%rho) >= c) value = transport%rho(c)
      end if
   end function lowmach_cell_rho_current


   real(rk) function lowmach_cell_rho_projection(fields, transport, c) result(value)
      type(flow_fields_t), intent(in) :: fields
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: c

      value = lowmach_cell_rho_current(transport, c)
      if (allocated(fields%projection_rho)) then
         if (size(fields%projection_rho) >= c) value = fields%projection_rho(c)
      end if
   end function lowmach_cell_rho_projection


   real(rk) function lowmach_cell_divu(mesh, fields, c) result(value)
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(in) :: fields
      integer, intent(in) :: c

      integer :: lf, fid
      real(rk) :: vol, flux_out

      value = 0.0_rk
      if (.not. allocated(fields%face_flux)) return
      if (c < 1 .or. c > mesh%ncells) return

      vol = mesh%cells(c)%volume
      if (vol <= 0.0_rk) return

      do lf = 1, mesh%ncell_faces(c)
         fid = mesh%cell_faces(lf, c)
         if (mesh%faces(fid)%owner == c) then
            flux_out = fields%face_flux(fid)
         else
            flux_out = -fields%face_flux(fid)
         end if
         value = value + flux_out / vol
      end do

   end function lowmach_cell_divu


   real(rk) function lowmach_cell_div_mass_flux(mesh, fields, c) result(value)
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(in) :: fields
      integer, intent(in) :: c

      integer :: lf, fid
      real(rk) :: vol, flux_out

      value = 0.0_rk
      if (.not. allocated(fields%mass_flux)) return
      if (c < 1 .or. c > mesh%ncells) return

      vol = mesh%cells(c)%volume
      if (vol <= 0.0_rk) return

      do lf = 1, mesh%ncell_faces(c)
         fid = mesh%cell_faces(lf, c)
         if (mesh%faces(fid)%owner == c) then
            flux_out = fields%mass_flux(fid)
         else
            flux_out = -fields%mass_flux(fid)
         end if
         value = value + flux_out / vol
      end do

   end function lowmach_cell_div_mass_flux


   real(rk) function lowmach_cell_u_dot_grad_rho(mesh, fields, transport, c) result(value)
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(in) :: fields
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: c

      integer :: lf, fid, nb
      real(rk) :: vol, flux_out, rho_c, rho_nb, rho_f

      value = 0.0_rk
      if (.not. allocated(fields%face_flux)) return
      if (.not. allocated(transport%rho)) return
      if (c < 1 .or. c > mesh%ncells) return
      if (size(transport%rho) < c) return

      vol = mesh%cells(c)%volume
      if (vol <= 0.0_rk) return

      rho_c = transport%rho(c)

      do lf = 1, mesh%ncell_faces(c)
         fid = mesh%cell_faces(lf, c)

         if (mesh%faces(fid)%owner == c) then
            flux_out = fields%face_flux(fid)
            nb = mesh%faces(fid)%neighbor
            if (nb <= 0) nb = mesh%faces(fid)%periodic_neighbor
         else
            flux_out = -fields%face_flux(fid)
            nb = mesh%faces(fid)%owner
         end if

         if (nb > 0 .and. nb <= mesh%ncells .and. size(transport%rho) >= nb) then
            rho_nb = transport%rho(nb)
            rho_f = 0.5_rk * (rho_c + rho_nb)
         else
            rho_f = rho_c
         end if

         value = value + flux_out * (rho_f - rho_c)
      end do

      value = value / vol

   end function lowmach_cell_u_dot_grad_rho

   !> Write aggregate species and enthalpy conservation trend diagnostics.
   !!
   !! This diagnostics-only routine is a first-pass conservation tracker for
   !! the variable-density path.  It reports global integrals,
   !! boundary-flux estimates, and step-to-step balance defects for:
   !!
   !!   - total mass
   !!   - transported species mass
   !!   - transported species sum
   !!   - rho*h
   !!
   !! Species boundary fluxes use stored face mass flux and owner-cell species
   !! values on physical boundary faces.  The rho*h flux is advective only; the
   !! aggregate enthalpy defect therefore excludes reconstructed conductive
   !! boundary fluxes and should be interpreted as a trend diagnostic, not a
   !! complete energy closure proof.
   subroutine write_species_energy_conservation_diagnostics(mesh, flow, params, fields, species, energy, transport, step, time)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(species_fields_t), intent(in) :: species
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: step
      real(rk), intent(in) :: time

      integer :: c, f, k, owner, nb, ierr, unit_id, nsp
      logical, save :: aggregate_initialized = .false.
      logical, save :: species_initialized = .false.
      logical, save :: previous_initialized = .false.
      character(len=1024) :: filename
      character(len=64) :: label
      real(rk) :: vol, rho_c, y_sum, h_c, t_c, mflux
      real(rk) :: dt_since_previous, denom
      real(rk) :: local_volume, volume_total
      real(rk) :: local_mass, total_mass
      real(rk) :: local_boundary_mass_flux, net_boundary_mass_flux
      real(rk) :: local_species_mass_sum, species_mass_sum
      real(rk) :: local_boundary_species_flux_sum, net_boundary_species_flux_sum
      real(rk) :: local_y_sum_min, y_sum_min
      real(rk) :: local_y_sum_max, y_sum_max
      real(rk) :: local_y_sum_int, y_sum_int
      real(rk) :: local_y_sum_l2_num, y_sum_l2_num
      real(rk) :: species_sum_mean, species_sum_l2
      real(rk) :: local_y_sum_minus_one_abs_max, y_sum_minus_one_abs_max
      real(rk) :: local_y_sum_minus_one_l2_num, y_sum_minus_one_l2_num
      real(rk) :: species_sum_minus_one_l2
      real(rk) :: local_rho_h_integral, rho_h_integral
      real(rk) :: local_boundary_rho_h_advective_flux, net_boundary_rho_h_advective_flux
      real(rk) :: local_qrad_integral, qrad_integral
      real(rk) :: local_rho_species_hdiff_integral, rho_species_hdiff_integral
      real(rk) :: local_h_min, h_min
      real(rk) :: local_h_max, h_max
      real(rk) :: local_t_min, t_min
      real(rk) :: local_t_max, t_max
      real(rk) :: delta_mass, mass_rate, mass_balance_defect, rel_mass_balance_defect
      real(rk) :: delta_rho_h, rho_h_rate, rho_h_advective_balance_defect_no_conduction
      real(rk), save :: previous_time = 0.0_rk
      real(rk), save :: previous_mass = 0.0_rk
      real(rk), save :: previous_rho_h = 0.0_rk
      real(rk), allocatable, save :: previous_species_mass(:)
      real(rk), allocatable :: local_species_mass(:), species_mass(:)
      real(rk), allocatable :: local_boundary_species_flux(:), boundary_species_flux(:)
      real(rk), allocatable :: species_delta_mass(:), species_mass_rate(:)
      real(rk), allocatable :: species_balance_defect(:), species_relative_balance_defect(:)
      real(rk) :: send_val, recv_val

      if (.not. params%write_diagnostics) return

      nsp = 0
      if (params%enable_species) nsp = max(0, species%nspecies)

      allocate(local_species_mass(max(1, nsp)))
      allocate(species_mass(max(1, nsp)))
      allocate(local_boundary_species_flux(max(1, nsp)))
      allocate(boundary_species_flux(max(1, nsp)))
      allocate(species_delta_mass(max(1, nsp)))
      allocate(species_mass_rate(max(1, nsp)))
      allocate(species_balance_defect(max(1, nsp)))
      allocate(species_relative_balance_defect(max(1, nsp)))

      local_species_mass = 0.0_rk
      species_mass = 0.0_rk
      local_boundary_species_flux = 0.0_rk
      boundary_species_flux = 0.0_rk
      species_delta_mass = 0.0_rk
      species_mass_rate = 0.0_rk
      species_balance_defect = 0.0_rk
      species_relative_balance_defect = 0.0_rk

      local_volume = 0.0_rk
      local_mass = 0.0_rk
      local_boundary_mass_flux = 0.0_rk
      local_species_mass_sum = 0.0_rk
      local_boundary_species_flux_sum = 0.0_rk
      local_y_sum_min = huge(1.0_rk)
      local_y_sum_max = -huge(1.0_rk)
      local_y_sum_int = 0.0_rk
      local_y_sum_l2_num = 0.0_rk
      local_y_sum_minus_one_abs_max = 0.0_rk
      local_y_sum_minus_one_l2_num = 0.0_rk
      local_rho_h_integral = 0.0_rk
      local_boundary_rho_h_advective_flux = 0.0_rk
      local_qrad_integral = 0.0_rk
      local_rho_species_hdiff_integral = 0.0_rk
      local_h_min = huge(1.0_rk)
      local_h_max = -huge(1.0_rk)
      local_t_min = huge(1.0_rk)
      local_t_max = -huge(1.0_rk)

      do c = 1, mesh%ncells
         if (.not. flow%owned(c)) cycle

         vol = mesh%cells(c)%volume
         if (vol <= 0.0_rk) cycle

         rho_c = params%rho
         if (allocated(transport%rho)) then
            if (size(transport%rho) >= c) rho_c = transport%rho(c)
         end if

         local_volume = local_volume + vol
         local_mass = local_mass + rho_c * vol

         y_sum = 0.0_rk
         if (nsp > 0 .and. allocated(species%Y)) then
            if (size(species%Y, 2) >= c) then
               do k = 1, nsp
                  if (size(species%Y, 1) >= k) then
                     local_species_mass(k) = local_species_mass(k) + rho_c * species%Y(k, c) * vol
                     y_sum = y_sum + species%Y(k, c)
                  end if
               end do
            end if
         end if

         if (nsp > 0) then
            local_species_mass_sum = local_species_mass_sum + rho_c * y_sum * vol
            local_y_sum_min = min(local_y_sum_min, y_sum)
            local_y_sum_max = max(local_y_sum_max, y_sum)
            local_y_sum_int = local_y_sum_int + y_sum * vol
            local_y_sum_l2_num = local_y_sum_l2_num + y_sum * y_sum * vol
            local_y_sum_minus_one_abs_max = max(local_y_sum_minus_one_abs_max, abs(y_sum - 1.0_rk))
            local_y_sum_minus_one_l2_num = local_y_sum_minus_one_l2_num + (y_sum - 1.0_rk)**2 * vol
         end if

         if (params%enable_energy .and. allocated(energy%h)) then
            if (size(energy%h) >= c) then
               h_c = energy%h(c)
               local_rho_h_integral = local_rho_h_integral + rho_c * h_c * vol
               local_h_min = min(local_h_min, h_c)
               local_h_max = max(local_h_max, h_c)
            end if
         end if

         if (params%enable_energy .and. allocated(energy%T)) then
            if (size(energy%T) >= c) then
               t_c = energy%T(c)
               local_t_min = min(local_t_min, t_c)
               local_t_max = max(local_t_max, t_c)
            end if
         end if

         if (params%enable_energy .and. allocated(energy%qrad)) then
            if (size(energy%qrad) >= c) local_qrad_integral = local_qrad_integral + energy%qrad(c) * vol
         end if

         if (params%enable_energy .and. allocated(energy%species_enthalpy_diffusion)) then
            if (size(energy%species_enthalpy_diffusion) >= c) then
               local_rho_species_hdiff_integral = local_rho_species_hdiff_integral + &
                  rho_c * energy%species_enthalpy_diffusion(c) * vol
            end if
         end if
      end do

      if (allocated(fields%mass_flux)) then
         do f = 1, mesh%nfaces
            owner = mesh%faces(f)%owner
            nb = mesh%faces(f)%neighbor
            if (nb <= 0) nb = mesh%faces(f)%periodic_neighbor

            if (nb <= 0 .and. owner >= 1 .and. owner <= mesh%ncells) then
               if (.not. flow%owned(owner)) cycle

               mflux = fields%mass_flux(f)
               local_boundary_mass_flux = local_boundary_mass_flux + mflux

               if (nsp > 0 .and. allocated(species%Y)) then
                  if (size(species%Y, 2) >= owner) then
                     do k = 1, nsp
                        if (size(species%Y, 1) >= k) then
                           local_boundary_species_flux(k) = local_boundary_species_flux(k) + mflux * species%Y(k, owner)
                           local_boundary_species_flux_sum = local_boundary_species_flux_sum + mflux * species%Y(k, owner)
                        end if
                     end do
                  end if
               end if

               if (params%enable_energy .and. allocated(energy%h)) then
                  if (size(energy%h) >= owner) then
                     local_boundary_rho_h_advective_flux = local_boundary_rho_h_advective_flux + mflux * energy%h(owner)
                  end if
               end if
            end if
         end do
      end if

      call reduce_sum(local_volume, volume_total, flow)
      call reduce_sum(local_mass, total_mass, flow)
      call reduce_sum(local_boundary_mass_flux, net_boundary_mass_flux, flow)
      call reduce_sum(local_species_mass_sum, species_mass_sum, flow)
      call reduce_sum(local_boundary_species_flux_sum, net_boundary_species_flux_sum, flow)
      call reduce_sum(local_y_sum_int, y_sum_int, flow)
      call reduce_sum(local_y_sum_l2_num, y_sum_l2_num, flow)
      call reduce_sum(local_y_sum_minus_one_l2_num, y_sum_minus_one_l2_num, flow)
      call reduce_sum(local_rho_h_integral, rho_h_integral, flow)
      call reduce_sum(local_boundary_rho_h_advective_flux, net_boundary_rho_h_advective_flux, flow)
      call reduce_sum(local_qrad_integral, qrad_integral, flow)
      call reduce_sum(local_rho_species_hdiff_integral, rho_species_hdiff_integral, flow)

      if (nsp > 0) then
         send_val = local_y_sum_min
         call MPI_Allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
         if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing species sum min')
         y_sum_min = recv_val

         send_val = local_y_sum_max
         call MPI_Allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
         if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing species sum max')
         y_sum_max = recv_val

         send_val = local_y_sum_minus_one_abs_max
         call MPI_Allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
         if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing species sum minus one max')
         y_sum_minus_one_abs_max = recv_val
      else
         y_sum_min = 0.0_rk
         y_sum_max = 0.0_rk
         y_sum_minus_one_abs_max = 0.0_rk
      end if

      if (params%enable_energy .and. allocated(energy%h)) then
         send_val = local_h_min
         call MPI_Allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
         if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing h min')
         h_min = recv_val

         send_val = local_h_max
         call MPI_Allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
         if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing h max')
         h_max = recv_val
      else
         h_min = 0.0_rk
         h_max = 0.0_rk
      end if

      if (params%enable_energy .and. allocated(energy%T)) then
         send_val = local_t_min
         call MPI_Allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
         if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing T min')
         t_min = recv_val

         send_val = local_t_max
         call MPI_Allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
         if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing T max')
         t_max = recv_val
      else
         t_min = 0.0_rk
         t_max = 0.0_rk
      end if

      if (nsp > 0) then
         call MPI_Allreduce(local_species_mass, species_mass, nsp, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
         if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing species mass integrals')

         call MPI_Allreduce(local_boundary_species_flux, boundary_species_flux, nsp, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
         if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing species boundary fluxes')
      end if

      if (volume_total > 0.0_rk) then
         species_sum_mean = y_sum_int / volume_total
         species_sum_l2 = sqrt(max(0.0_rk, y_sum_l2_num / volume_total))
         species_sum_minus_one_l2 = sqrt(max(0.0_rk, y_sum_minus_one_l2_num / volume_total))
      else
         species_sum_mean = 0.0_rk
         species_sum_l2 = 0.0_rk
         species_sum_minus_one_l2 = 0.0_rk
      end if

      if (previous_initialized) then
         dt_since_previous = time - previous_time
      else
         dt_since_previous = 0.0_rk
      end if

      if (previous_initialized .and. dt_since_previous > tiny(1.0_rk)) then
         delta_mass = total_mass - previous_mass
         mass_rate = delta_mass / dt_since_previous
         mass_balance_defect = mass_rate + net_boundary_mass_flux
         denom = max(abs(net_boundary_mass_flux), abs(mass_rate), tiny(1.0_rk))
         rel_mass_balance_defect = abs(mass_balance_defect) / denom

         delta_rho_h = rho_h_integral - previous_rho_h
         rho_h_rate = delta_rho_h / dt_since_previous
         rho_h_advective_balance_defect_no_conduction = rho_h_rate + net_boundary_rho_h_advective_flux - &
                                                        qrad_integral - rho_species_hdiff_integral

         if (nsp > 0) then
            if (.not. allocated(previous_species_mass)) allocate(previous_species_mass(nsp))
            if (size(previous_species_mass) /= nsp) then
               deallocate(previous_species_mass)
               allocate(previous_species_mass(nsp))
               previous_species_mass = species_mass
            end if

            do k = 1, nsp
               species_delta_mass(k) = species_mass(k) - previous_species_mass(k)
               species_mass_rate(k) = species_delta_mass(k) / dt_since_previous
               species_balance_defect(k) = species_mass_rate(k) + boundary_species_flux(k)
               denom = max(abs(species_mass_rate(k)), abs(boundary_species_flux(k)), tiny(1.0_rk))
               species_relative_balance_defect(k) = abs(species_balance_defect(k)) / denom
            end do
         end if
      else
         delta_mass = 0.0_rk
         mass_rate = 0.0_rk
         mass_balance_defect = 0.0_rk
         rel_mass_balance_defect = 0.0_rk
         delta_rho_h = 0.0_rk
         rho_h_rate = 0.0_rk
         rho_h_advective_balance_defect_no_conduction = 0.0_rk
      end if

      if (flow%rank == 0) then
         call execute_command_line('mkdir -p ' // trim(params%output_dir) // '/diagnostics')

         filename = trim(params%output_dir) // '/diagnostics/species_energy_conservation.csv'
         if (.not. aggregate_initialized) then
            open(newunit=unit_id, file=trim(filename), status='replace', action='write')
            aggregate_initialized = .true.
            write(unit_id,'(a)') 'step,time,volume_total,total_mass,net_boundary_mass_flux,' // &
               'delta_time_since_previous,delta_mass_since_previous,mass_rate_since_previous,' // &
               'mass_balance_defect_since_previous,relative_mass_balance_defect_since_previous,' // &
               'transported_species_mass_sum,net_boundary_species_mass_flux_sum,' // &
               'species_sum_min,species_sum_max,species_sum_mean,species_sum_l2,' // &
               'species_sum_minus_one_abs_max,species_sum_minus_one_l2,' // &
               'rho_h_integral,net_boundary_rho_h_advective_flux,' // &
               'qrad_integral,rho_species_enthalpy_diffusion_integral,' // &
               'delta_rho_h_since_previous,rho_h_rate_since_previous,' // &
               'rho_h_advective_balance_defect_no_conduction,h_min,h_max,T_min,T_max'
         else
            open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
         end if

         call write_csv_value_i(unit_id, step)
         call write_csv_value_r(unit_id, time)
         call write_csv_value_r(unit_id, volume_total)
         call write_csv_value_r(unit_id, total_mass)
         call write_csv_value_r(unit_id, net_boundary_mass_flux)
         call write_csv_value_r(unit_id, dt_since_previous)
         call write_csv_value_r(unit_id, delta_mass)
         call write_csv_value_r(unit_id, mass_rate)
         call write_csv_value_r(unit_id, mass_balance_defect)
         call write_csv_value_r(unit_id, rel_mass_balance_defect)
         call write_csv_value_r(unit_id, species_mass_sum)
         call write_csv_value_r(unit_id, net_boundary_species_flux_sum)
         call write_csv_value_r(unit_id, y_sum_min)
         call write_csv_value_r(unit_id, y_sum_max)
         call write_csv_value_r(unit_id, species_sum_mean)
         call write_csv_value_r(unit_id, species_sum_l2)
         call write_csv_value_r(unit_id, y_sum_minus_one_abs_max)
         call write_csv_value_r(unit_id, species_sum_minus_one_l2)
         call write_csv_value_r(unit_id, rho_h_integral)
         call write_csv_value_r(unit_id, net_boundary_rho_h_advective_flux)
         call write_csv_value_r(unit_id, qrad_integral)
         call write_csv_value_r(unit_id, rho_species_hdiff_integral)
         call write_csv_value_r(unit_id, delta_rho_h)
         call write_csv_value_r(unit_id, rho_h_rate)
         call write_csv_value_r(unit_id, rho_h_advective_balance_defect_no_conduction)
         call write_csv_value_r(unit_id, h_min)
         call write_csv_value_r(unit_id, h_max)
         call write_csv_value_r(unit_id, t_min)
         call write_csv_value_r(unit_id, t_max)
         write(unit_id,*)
         close(unit_id)

         filename = trim(params%output_dir) // '/diagnostics/species_integrals.csv'
         if (nsp > 0) then
            if (.not. species_initialized) then
               open(newunit=unit_id, file=trim(filename), status='replace', action='write')
               species_initialized = .true.
               write(unit_id,'(a)', advance='no') 'step,time'
               do k = 1, nsp
                  write(label,'("Y",i0)') k
                  call write_csv_header(unit_id, trim(label)//'_mass_integral')
                  call write_csv_header(unit_id, trim(label)//'_net_boundary_mass_flux')
                  call write_csv_header(unit_id, trim(label)//'_delta_mass_since_previous')
                  call write_csv_header(unit_id, trim(label)//'_mass_rate_since_previous')
                  call write_csv_header(unit_id, trim(label)//'_balance_defect_since_previous')
                  call write_csv_header(unit_id, trim(label)//'_relative_balance_defect_since_previous')
               end do
               write(unit_id,*)
            else
               open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
            end if

            call write_csv_value_i(unit_id, step)
            call write_csv_value_r(unit_id, time)
            do k = 1, nsp
               call write_csv_value_r(unit_id, species_mass(k))
               call write_csv_value_r(unit_id, boundary_species_flux(k))
               call write_csv_value_r(unit_id, species_delta_mass(k))
               call write_csv_value_r(unit_id, species_mass_rate(k))
               call write_csv_value_r(unit_id, species_balance_defect(k))
               call write_csv_value_r(unit_id, species_relative_balance_defect(k))
            end do
            write(unit_id,*)
            close(unit_id)
         end if
      end if

      previous_time = time
      previous_mass = total_mass
      previous_rho_h = rho_h_integral
      previous_initialized = .true.

      if (nsp > 0) then
         if (.not. allocated(previous_species_mass)) allocate(previous_species_mass(nsp))
         if (size(previous_species_mass) /= nsp) then
            deallocate(previous_species_mass)
            allocate(previous_species_mass(nsp))
         end if
         previous_species_mass = species_mass(1:nsp)
      end if

      deallocate(local_species_mass)
      deallocate(species_mass)
      deallocate(local_boundary_species_flux)
      deallocate(boundary_species_flux)
      deallocate(species_delta_mass)
      deallocate(species_mass_rate)
      deallocate(species_balance_defect)
      deallocate(species_relative_balance_defect)

   end subroutine write_species_energy_conservation_diagnostics


   subroutine reduce_sum(local_value, global_value, flow)
      real(rk), intent(in) :: local_value
      real(rk), intent(out) :: global_value
      type(flow_mpi_t), intent(in) :: flow

      integer :: ierr
      real(rk) :: send_val, recv_val

      send_val = local_value
      call MPI_Allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure in reduce_sum')
      global_value = recv_val

   end subroutine reduce_sum


   subroutine write_csv_value_i(unit_id, value)
      integer, intent(in) :: unit_id
      integer, intent(in) :: value

      write(unit_id,'(i0)', advance='no') value

   end subroutine write_csv_value_i


   subroutine write_csv_value_r(unit_id, value)
      integer, intent(in) :: unit_id
      real(rk), intent(in) :: value

      write(unit_id,'(",",es16.8)', advance='no') value

   end subroutine write_csv_value_r


   subroutine write_csv_header(unit_id, text)
      integer, intent(in) :: unit_id
      character(len=*), intent(in) :: text

      write(unit_id,'(",",a)', advance='no') trim(text)

   end subroutine write_csv_header

   !> Write dedicated enthalpy/energy budget diagnostics.
   !!
   !! This diagnostics-only routine separates the terms that are currently
   !! available from the missing conductive-boundary-flux closure term.
   !!
   !! Sign convention:
   !!   boundary fluxes are positive outward from the domain.
   !!
   !! Reported balance without conduction:
   !!   d/dt integral(rho h dV)
   !! + net_boundary_advective_rho_h_flux
   !! - qrad_integral
   !! - rho_species_enthalpy_diffusion_integral
   !!
   !! If the full enthalpy equation is written with an outward conductive flux
   !! on the left-hand side, the conductive boundary flux required for closure is:
   !!   required_conductive_out = -balance_without_conduction
   subroutine write_enthalpy_energy_budget_diagnostics(mesh, flow, params, fields, energy, transport, step, time)
      type(mesh_t), intent(in) :: mesh
      type(flow_mpi_t), intent(in) :: flow
      type(case_params_t), intent(in) :: params
      type(flow_fields_t), intent(in) :: fields
      type(energy_fields_t), intent(in) :: energy
      type(transport_properties_t), intent(in) :: transport
      integer, intent(in) :: step
      real(rk), intent(in) :: time

      integer :: c, f, owner, nb, ierr, unit_id
      logical, save :: enthalpy_budget_initialized = .false.
      logical, save :: previous_initialized = .false.
      character(len=1024) :: filename
      real(rk), save :: previous_time = 0.0_rk
      real(rk), save :: previous_rho_h_integral = 0.0_rk
      real(rk) :: vol, rho_c, h_c, T_c, cp_c, lambda_c, mflux
      real(rk) :: dt_since_previous, delta_rho_h, rho_h_rate
      real(rk) :: balance_without_conduction, balance_with_reported_conduction
      real(rk) :: required_conductive_boundary_flux_out
      real(rk) :: relative_balance_without_conduction
      real(rk) :: denom
      real(rk) :: local_volume, volume_total
      real(rk) :: local_rho_h_integral, rho_h_integral
      real(rk) :: local_boundary_rho_h_advective_flux, net_boundary_rho_h_advective_flux
      real(rk) :: local_boundary_rho_h_advective_outflow, boundary_rho_h_advective_outflow
      real(rk) :: local_boundary_rho_h_advective_inflow, boundary_rho_h_advective_inflow
      real(rk) :: local_qrad_integral, qrad_integral
      real(rk) :: local_rho_species_hdiff_integral, rho_species_hdiff_integral
      real(rk) :: local_h_min, h_min
      real(rk) :: local_h_max, h_max
      real(rk) :: local_T_min, T_min
      real(rk) :: local_T_max, T_max
      real(rk) :: local_rho_min, rho_min
      real(rk) :: local_rho_max, rho_max
      real(rk) :: local_cp_min, cp_min
      real(rk) :: local_cp_max, cp_max
      real(rk) :: local_lambda_min, lambda_min
      real(rk) :: local_lambda_max, lambda_max
      real(rk) :: local_h_int, h_mean
      real(rk) :: local_T_int, T_mean
      real(rk) :: local_rho_int, rho_mean
      real(rk) :: local_cp_int, cp_mean
      real(rk) :: local_lambda_int, lambda_mean
      real(rk) :: local_rho_h_l2_num, rho_h_l2
      real(rk) :: sum_send(24), sum_recv(24)
      real(rk) :: min_send(5), min_recv(5)
      real(rk) :: max_send(5), max_recv(5)
      real(rk) :: reported_conductive_boundary_flux_out
      integer :: conductive_boundary_flux_available
      real(rk), save :: previous_cumulative_advective_boundary_rho_h_flux = 0.0_rk
      real(rk), save :: previous_cumulative_conductive_boundary_flux = 0.0_rk
      real(rk), save :: previous_cumulative_qrad_integral = 0.0_rk
      real(rk), save :: previous_cumulative_rho_species_hdiff_integral = 0.0_rk
      real(rk) :: cumulative_advective_boundary_rho_h_flux
      real(rk) :: cumulative_conductive_boundary_flux
      real(rk) :: cumulative_qrad_integral_total
      real(rk) :: cumulative_rho_species_hdiff_integral_total
      integer :: cumulative_energy_budget_available
      real(rk) :: cumulative_advective_boundary_rho_h_flux_average
      real(rk) :: cumulative_conductive_boundary_flux_average
      real(rk) :: cumulative_qrad_integral_average
      real(rk) :: cumulative_rho_species_hdiff_integral_average
      real(rk) :: cumulative_operator_balance_defect
      real(rk) :: relative_cumulative_operator_balance_defect
      real(rk), save :: previous_operator_consistent_rho_h_integral = 0.0_rk
      real(rk) :: operator_consistent_rho_h_integral
      real(rk) :: operator_consistent_delta_rho_h
      real(rk) :: operator_consistent_rho_h_rate
      real(rk) :: operator_consistent_cumulative_balance_defect
      real(rk) :: relative_operator_consistent_cumulative_balance_defect
      real(rk), save :: previous_cumulative_energy_update_delta_integral = 0.0_rk
      real(rk), save :: previous_cumulative_energy_update_rhs_integral = 0.0_rk
      real(rk) :: cumulative_energy_update_delta_rate_average
      real(rk) :: cumulative_energy_update_rhs_rate_average
      real(rk) :: output_state_density_reconciliation_rate
      real(rk) :: operator_consistent_density_reconciliation_rate
      real(rk) :: output_state_budget_defect_after_density_reconciliation
      real(rk) :: operator_consistent_budget_defect_after_density_reconciliation
      real(rk) :: rel_output_recon_defect
      real(rk) :: rel_operator_recon_defect

      if (.not. params%write_diagnostics) return
      if (.not. params%enable_energy) return
      if (.not. allocated(energy%h)) return
      if (.not. allocated(transport%rho)) return

      local_volume = 0.0_rk
      local_rho_h_integral = 0.0_rk
      local_boundary_rho_h_advective_flux = 0.0_rk
      local_boundary_rho_h_advective_outflow = 0.0_rk
      local_boundary_rho_h_advective_inflow = 0.0_rk
      local_qrad_integral = 0.0_rk
      local_rho_species_hdiff_integral = 0.0_rk
      local_h_min = huge(1.0_rk)
      local_h_max = -huge(1.0_rk)
      local_T_min = huge(1.0_rk)
      local_T_max = -huge(1.0_rk)
      local_rho_min = huge(1.0_rk)
      local_rho_max = -huge(1.0_rk)
      local_cp_min = huge(1.0_rk)
      local_cp_max = -huge(1.0_rk)
      local_lambda_min = huge(1.0_rk)
      local_lambda_max = -huge(1.0_rk)
      local_h_int = 0.0_rk
      local_T_int = 0.0_rk
      local_rho_int = 0.0_rk
      local_cp_int = 0.0_rk
      local_lambda_int = 0.0_rk
      local_rho_h_l2_num = 0.0_rk

      do c = 1, mesh%ncells
         if (.not. flow%owned(c)) cycle

         vol = mesh%cells(c)%volume
         if (vol <= 0.0_rk) cycle

         rho_c = transport%rho(c)
         h_c = energy%h(c)

         T_c = 0.0_rk
         if (allocated(energy%T)) then
            if (size(energy%T) >= c) T_c = energy%T(c)
         end if

         cp_c = 0.0_rk
         if (allocated(energy%cp)) then
            if (size(energy%cp) >= c) cp_c = energy%cp(c)
         end if

         lambda_c = 0.0_rk
         if (allocated(energy%lambda)) then
            if (size(energy%lambda) >= c) lambda_c = energy%lambda(c)
         end if

         local_volume = local_volume + vol
         local_rho_h_integral = local_rho_h_integral + rho_c * h_c * vol
         local_rho_h_l2_num = local_rho_h_l2_num + (rho_c * h_c) * (rho_c * h_c) * vol

         local_h_min = min(local_h_min, h_c)
         local_h_max = max(local_h_max, h_c)
         local_T_min = min(local_T_min, T_c)
         local_T_max = max(local_T_max, T_c)
         local_rho_min = min(local_rho_min, rho_c)
         local_rho_max = max(local_rho_max, rho_c)
         local_cp_min = min(local_cp_min, cp_c)
         local_cp_max = max(local_cp_max, cp_c)
         local_lambda_min = min(local_lambda_min, lambda_c)
         local_lambda_max = max(local_lambda_max, lambda_c)

         local_h_int = local_h_int + h_c * vol
         local_T_int = local_T_int + T_c * vol
         local_rho_int = local_rho_int + rho_c * vol
         local_cp_int = local_cp_int + cp_c * vol
         local_lambda_int = local_lambda_int + lambda_c * vol

         if (allocated(energy%qrad)) then
            if (size(energy%qrad) >= c) local_qrad_integral = local_qrad_integral + energy%qrad(c) * vol
         end if

         if (allocated(energy%species_enthalpy_diffusion)) then
            if (size(energy%species_enthalpy_diffusion) >= c) then
               local_rho_species_hdiff_integral = local_rho_species_hdiff_integral + &
                  rho_c * energy%species_enthalpy_diffusion(c) * vol
            end if
         end if
      end do

      if (allocated(fields%mass_flux)) then
         do f = 1, mesh%nfaces
            owner = mesh%faces(f)%owner
            nb = mesh%faces(f)%neighbor
            if (nb <= 0) nb = mesh%faces(f)%periodic_neighbor

            if (nb <= 0 .and. owner >= 1 .and. owner <= mesh%ncells) then
               if (.not. flow%owned(owner)) cycle
               if (size(energy%h) < owner) cycle

               mflux = fields%mass_flux(f)
               h_c = energy%h(owner)

               local_boundary_rho_h_advective_flux = local_boundary_rho_h_advective_flux + mflux * h_c
               if (mflux >= 0.0_rk) then
                  local_boundary_rho_h_advective_outflow = local_boundary_rho_h_advective_outflow + mflux * h_c
               else
                  local_boundary_rho_h_advective_inflow = local_boundary_rho_h_advective_inflow + mflux * h_c
               end if
            end if
         end do
      end if

      sum_send = 0.0_rk
      sum_send(1) = local_volume
      sum_send(2) = local_rho_h_integral
      sum_send(3) = local_boundary_rho_h_advective_flux
      sum_send(4) = local_boundary_rho_h_advective_outflow
      sum_send(5) = local_boundary_rho_h_advective_inflow
      sum_send(6) = local_qrad_integral
      sum_send(7) = local_rho_species_hdiff_integral
      sum_send(8) = local_h_int
      sum_send(9) = local_T_int
      sum_send(10) = local_rho_int
      sum_send(11) = local_cp_int
      sum_send(12) = local_lambda_int
      sum_send(13) = local_rho_h_l2_num
      sum_send(14) = energy%last_conductive_boundary_flux_out
      sum_send(15) = real(energy%last_conductive_boundary_flux_available, rk)
      sum_send(16) = energy%cumulative_boundary_rho_h_advective_flux_out
      sum_send(17) = energy%cumulative_boundary_rho_h_conductive_flux_out
      sum_send(18) = energy%cumulative_qrad_integral
      sum_send(19) = energy%cumulative_rho_species_hdiff_integral
      sum_send(20) = real(energy%cumulative_energy_budget_available, rk)
      call MPI_Allreduce(sum_send, sum_recv, 24, MPI_DOUBLE_PRECISION, MPI_SUM, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing enthalpy budget sums')

      min_send(1) = local_h_min
      min_send(2) = local_T_min
      min_send(3) = local_rho_min
      min_send(4) = local_cp_min
      min_send(5) = local_lambda_min
      call MPI_Allreduce(min_send, min_recv, 5, MPI_DOUBLE_PRECISION, MPI_MIN, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing enthalpy budget minima')

      max_send(1) = local_h_max
      max_send(2) = local_T_max
      max_send(3) = local_rho_max
      max_send(4) = local_cp_max
      max_send(5) = local_lambda_max
      call MPI_Allreduce(max_send, max_recv, 5, MPI_DOUBLE_PRECISION, MPI_MAX, flow%comm, ierr)
      if (ierr /= MPI_SUCCESS) call fatal_error('output', 'MPI failure reducing enthalpy budget maxima')

      volume_total = sum_recv(1)
      rho_h_integral = sum_recv(2)
      net_boundary_rho_h_advective_flux = sum_recv(3)
      boundary_rho_h_advective_outflow = sum_recv(4)
      boundary_rho_h_advective_inflow = sum_recv(5)
      qrad_integral = sum_recv(6)
      rho_species_hdiff_integral = sum_recv(7)

      if (volume_total > 0.0_rk) then
         h_mean = sum_recv(8) / volume_total
         T_mean = sum_recv(9) / volume_total
         rho_mean = sum_recv(10) / volume_total
         cp_mean = sum_recv(11) / volume_total
         lambda_mean = sum_recv(12) / volume_total
         rho_h_l2 = sqrt(max(0.0_rk, sum_recv(13) / volume_total))
      else
         h_mean = 0.0_rk
         T_mean = 0.0_rk
         rho_mean = 0.0_rk
         cp_mean = 0.0_rk
         lambda_mean = 0.0_rk
         rho_h_l2 = 0.0_rk
      end if

      h_min = min_recv(1)
      T_min = min_recv(2)
      rho_min = min_recv(3)
      cp_min = min_recv(4)
      lambda_min = min_recv(5)
      h_max = max_recv(1)
      T_max = max_recv(2)
      rho_max = max_recv(3)
      cp_max = max_recv(4)
      lambda_max = max_recv(5)

      if (previous_initialized) then
         dt_since_previous = time - previous_time
      else
         dt_since_previous = 0.0_rk
      end if

      if (previous_initialized .and. dt_since_previous > tiny(1.0_rk)) then
         delta_rho_h = rho_h_integral - previous_rho_h_integral
         rho_h_rate = delta_rho_h / dt_since_previous
      else
         delta_rho_h = 0.0_rk
         rho_h_rate = 0.0_rk
      end if

      reported_conductive_boundary_flux_out = sum_recv(14)
      if (sum_recv(15) > 0.5_rk) then
         conductive_boundary_flux_available = 1
      else
         conductive_boundary_flux_available = 0
      end if

      cumulative_advective_boundary_rho_h_flux = sum_recv(16)
      cumulative_conductive_boundary_flux = sum_recv(17)
      cumulative_qrad_integral_total = sum_recv(18)
      cumulative_rho_species_hdiff_integral_total = sum_recv(19)
      if (sum_recv(20) > 0.5_rk) then
         cumulative_energy_budget_available = 1
      else
         cumulative_energy_budget_available = 0
      end if

            if (previous_initialized .and. dt_since_previous > tiny(1.0_rk) .and. &
          cumulative_energy_budget_available == 1) then
         cumulative_advective_boundary_rho_h_flux_average = &
            (cumulative_advective_boundary_rho_h_flux - previous_cumulative_advective_boundary_rho_h_flux) / &
            dt_since_previous
         cumulative_conductive_boundary_flux_average = &
            (cumulative_conductive_boundary_flux - previous_cumulative_conductive_boundary_flux) / &
            dt_since_previous
         cumulative_qrad_integral_average = &
            (cumulative_qrad_integral_total - previous_cumulative_qrad_integral) / dt_since_previous
         cumulative_rho_species_hdiff_integral_average = &
            (cumulative_rho_species_hdiff_integral_total - previous_cumulative_rho_species_hdiff_integral) / &
            dt_since_previous
         cumulative_operator_balance_defect = rho_h_rate + &
            cumulative_advective_boundary_rho_h_flux_average + &
            cumulative_conductive_boundary_flux_average - &
            cumulative_qrad_integral_average - &
            cumulative_rho_species_hdiff_integral_average
         denom = max(abs(rho_h_rate), abs(cumulative_advective_boundary_rho_h_flux_average), &
                     abs(cumulative_conductive_boundary_flux_average), abs(cumulative_qrad_integral_average), &
                     abs(cumulative_rho_species_hdiff_integral_average), tiny(1.0_rk))
         relative_cumulative_operator_balance_defect = abs(cumulative_operator_balance_defect) / denom
      else
         cumulative_advective_boundary_rho_h_flux_average = 0.0_rk
         cumulative_conductive_boundary_flux_average = 0.0_rk
         cumulative_qrad_integral_average = 0.0_rk
         cumulative_rho_species_hdiff_integral_average = 0.0_rk
         cumulative_operator_balance_defect = 0.0_rk
         relative_cumulative_operator_balance_defect = 0.0_rk
      end if

            operator_consistent_rho_h_integral = energy%last_operator_consistent_rho_h_integral
      if (previous_initialized .and. dt_since_previous > tiny(1.0_rk) .and. &
          cumulative_energy_budget_available == 1) then
         operator_consistent_delta_rho_h = operator_consistent_rho_h_integral - &
                                           previous_operator_consistent_rho_h_integral
         operator_consistent_rho_h_rate = operator_consistent_delta_rho_h / dt_since_previous
         operator_consistent_cumulative_balance_defect = operator_consistent_rho_h_rate + &
            cumulative_advective_boundary_rho_h_flux_average + &
            cumulative_conductive_boundary_flux_average - &
            cumulative_qrad_integral_average - &
            cumulative_rho_species_hdiff_integral_average
         denom = max(abs(operator_consistent_rho_h_rate), &
                     abs(cumulative_advective_boundary_rho_h_flux_average), &
                     abs(cumulative_conductive_boundary_flux_average), &
                     abs(cumulative_qrad_integral_average), &
                     abs(cumulative_rho_species_hdiff_integral_average), tiny(1.0_rk))
         relative_operator_consistent_cumulative_balance_defect = &
            abs(operator_consistent_cumulative_balance_defect) / denom
      else
         operator_consistent_delta_rho_h = 0.0_rk
         operator_consistent_rho_h_rate = 0.0_rk
         operator_consistent_cumulative_balance_defect = 0.0_rk
         relative_operator_consistent_cumulative_balance_defect = 0.0_rk
      end if

            if (previous_initialized .and. dt_since_previous > tiny(1.0_rk) .and. &
          cumulative_energy_budget_available == 1) then
         cumulative_energy_update_delta_rate_average = &
            (energy%cumulative_energy_update_delta_integral - previous_cumulative_energy_update_delta_integral) / &
            dt_since_previous
         cumulative_energy_update_rhs_rate_average = &
            (energy%cumulative_energy_update_rhs_integral - previous_cumulative_energy_update_rhs_integral) / &
            dt_since_previous
         output_state_density_reconciliation_rate = rho_h_rate - cumulative_energy_update_delta_rate_average
         operator_consistent_density_reconciliation_rate = operator_consistent_rho_h_rate - &
                                                     cumulative_energy_update_delta_rate_average
         output_state_budget_defect_after_density_reconciliation = &
            cumulative_operator_balance_defect - output_state_density_reconciliation_rate
         operator_consistent_budget_defect_after_density_reconciliation = &
            operator_consistent_cumulative_balance_defect - operator_consistent_density_reconciliation_rate
         denom = max(abs(cumulative_energy_update_delta_rate_average), &
                     abs(cumulative_energy_update_rhs_rate_average), &
                     abs(cumulative_advective_boundary_rho_h_flux_average), &
                     abs(cumulative_conductive_boundary_flux_average), &
                     abs(cumulative_qrad_integral_average), &
                     abs(cumulative_rho_species_hdiff_integral_average), tiny(1.0_rk))
         rel_output_recon_defect = &
            abs(output_state_budget_defect_after_density_reconciliation) / denom
         rel_operator_recon_defect = &
            abs(operator_consistent_budget_defect_after_density_reconciliation) / denom
      else
         cumulative_energy_update_delta_rate_average = 0.0_rk
         cumulative_energy_update_rhs_rate_average = 0.0_rk
         output_state_density_reconciliation_rate = 0.0_rk
         operator_consistent_density_reconciliation_rate = 0.0_rk
         output_state_budget_defect_after_density_reconciliation = 0.0_rk
         operator_consistent_budget_defect_after_density_reconciliation = 0.0_rk
         rel_output_recon_defect = 0.0_rk
         rel_operator_recon_defect = 0.0_rk
      end if

      balance_without_conduction = rho_h_rate + net_boundary_rho_h_advective_flux - &
                                    qrad_integral - rho_species_hdiff_integral
      required_conductive_boundary_flux_out = -balance_without_conduction
      balance_with_reported_conduction = balance_without_conduction + reported_conductive_boundary_flux_out

      denom = max(abs(rho_h_rate), abs(net_boundary_rho_h_advective_flux), abs(qrad_integral), &
                  abs(rho_species_hdiff_integral), tiny(1.0_rk))
      relative_balance_without_conduction = abs(balance_without_conduction) / denom

      if (flow%rank == 0) then
         call execute_command_line('mkdir -p ' // trim(params%output_dir) // '/diagnostics')

         filename = trim(params%output_dir) // '/diagnostics/enthalpy_energy_budget.csv'
         if (.not. enthalpy_budget_initialized) then
            open(newunit=unit_id, file=trim(filename), status='replace', action='write')
            enthalpy_budget_initialized = .true.
            write(unit_id,'(a)') 'step,time,volume_total,rho_h_integral,delta_time_since_previous,' // &
               'delta_rho_h_since_previous,rho_h_rate_since_previous,' // &
               'net_boundary_rho_h_advective_flux,boundary_rho_h_advective_outflow,' // &
               'boundary_rho_h_advective_inflow,qrad_integral,' // &
               'rho_species_enthalpy_diffusion_integral,reported_conductive_boundary_flux_out,' // &
               'conductive_boundary_flux_available,required_conductive_boundary_flux_out,' // &
               'balance_defect_without_conduction,balance_defect_with_reported_conduction,' // &
               'relative_balance_defect_without_conduction,' // &
               'cumulative_advective_boundary_rho_h_flux_average,' // &
               'cumulative_conductive_boundary_flux_average,cumulative_qrad_integral_average,' // &
               'cumulative_rho_species_enthalpy_diffusion_integral_average,' // &
               'cumulative_energy_budget_available,cumulative_operator_balance_defect,' // &
               'relative_cumulative_operator_balance_defect,' // &
               'last_energy_update_delta_rate_integral,last_energy_update_rhs_integral,' // &
               'last_energy_update_balance_defect,relative_last_energy_update_balance_defect,' // &
               'operator_consistent_rho_h_integral,' // &
               'operator_consistent_delta_rho_h_since_previous,' // &
               'operator_consistent_rho_h_rate_since_previous,' // &
               'operator_consistent_cumulative_balance_defect,' // &
               'relative_operator_consistent_cumulative_balance_defect,' // &
               'cumulative_energy_update_delta_rate_average,' // &
               'cumulative_energy_update_rhs_rate_average,' // &
               'output_state_density_reconciliation_rate,' // &
               'operator_consistent_density_reconciliation_rate,' // &
               'output_state_budget_defect_after_density_reconciliation,' // &
               'operator_consistent_budget_defect_after_density_reconciliation,' // &
               'rel_output_recon_defect,' // &
               'rel_operator_recon_defect,' // &
               'h_min,h_max,h_mean,T_min,T_max,T_mean,' // &
               'rho_min,rho_max,rho_mean,cp_min,cp_max,cp_mean,lambda_min,lambda_max,lambda_mean,rho_h_l2'
         else
            open(newunit=unit_id, file=trim(filename), status='unknown', position='append', action='write')
         end if

         call write_csv_value_i(unit_id, step)
         call write_csv_value_r(unit_id, time)
         call write_csv_value_r(unit_id, volume_total)
         call write_csv_value_r(unit_id, rho_h_integral)
         call write_csv_value_r(unit_id, dt_since_previous)
         call write_csv_value_r(unit_id, delta_rho_h)
         call write_csv_value_r(unit_id, rho_h_rate)
         call write_csv_value_r(unit_id, net_boundary_rho_h_advective_flux)
         call write_csv_value_r(unit_id, boundary_rho_h_advective_outflow)
         call write_csv_value_r(unit_id, boundary_rho_h_advective_inflow)
         call write_csv_value_r(unit_id, qrad_integral)
         call write_csv_value_r(unit_id, rho_species_hdiff_integral)
         call write_csv_value_r(unit_id, reported_conductive_boundary_flux_out)
         write(unit_id,'(",",i0)', advance='no') conductive_boundary_flux_available
         call write_csv_value_r(unit_id, required_conductive_boundary_flux_out)
         call write_csv_value_r(unit_id, balance_without_conduction)
         call write_csv_value_r(unit_id, balance_with_reported_conduction)
         call write_csv_value_r(unit_id, relative_balance_without_conduction)
         call write_csv_value_r(unit_id, cumulative_advective_boundary_rho_h_flux_average)
         call write_csv_value_r(unit_id, cumulative_conductive_boundary_flux_average)
         call write_csv_value_r(unit_id, cumulative_qrad_integral_average)
         call write_csv_value_r(unit_id, cumulative_rho_species_hdiff_integral_average)
         write(unit_id,'(",",i0)', advance='no') cumulative_energy_budget_available
         call write_csv_value_r(unit_id, cumulative_operator_balance_defect)
         call write_csv_value_r(unit_id, relative_cumulative_operator_balance_defect)
         call write_csv_value_r(unit_id, energy%last_energy_update_delta_rate_integral)
         call write_csv_value_r(unit_id, energy%last_energy_update_rhs_integral)
         call write_csv_value_r(unit_id, energy%last_energy_update_balance_defect)
         call write_csv_value_r(unit_id, energy%relative_last_energy_update_balance_defect)
         call write_csv_value_r(unit_id, operator_consistent_rho_h_integral)
         call write_csv_value_r(unit_id, operator_consistent_delta_rho_h)
         call write_csv_value_r(unit_id, operator_consistent_rho_h_rate)
         call write_csv_value_r(unit_id, operator_consistent_cumulative_balance_defect)
         call write_csv_value_r(unit_id, relative_operator_consistent_cumulative_balance_defect)
         call write_csv_value_r(unit_id, cumulative_energy_update_delta_rate_average)
         call write_csv_value_r(unit_id, cumulative_energy_update_rhs_rate_average)
         call write_csv_value_r(unit_id, output_state_density_reconciliation_rate)
         call write_csv_value_r(unit_id, operator_consistent_density_reconciliation_rate)
         call write_csv_value_r(unit_id, output_state_budget_defect_after_density_reconciliation)
         call write_csv_value_r(unit_id, operator_consistent_budget_defect_after_density_reconciliation)
         call write_csv_value_r(unit_id, rel_output_recon_defect)
         call write_csv_value_r(unit_id, rel_operator_recon_defect)
         call write_csv_value_r(unit_id, h_min)
         call write_csv_value_r(unit_id, h_max)
         call write_csv_value_r(unit_id, h_mean)
         call write_csv_value_r(unit_id, T_min)
         call write_csv_value_r(unit_id, T_max)
         call write_csv_value_r(unit_id, T_mean)
         call write_csv_value_r(unit_id, rho_min)
         call write_csv_value_r(unit_id, rho_max)
         call write_csv_value_r(unit_id, rho_mean)
         call write_csv_value_r(unit_id, cp_min)
         call write_csv_value_r(unit_id, cp_max)
         call write_csv_value_r(unit_id, cp_mean)
         call write_csv_value_r(unit_id, lambda_min)
         call write_csv_value_r(unit_id, lambda_max)
         call write_csv_value_r(unit_id, lambda_mean)
         call write_csv_value_r(unit_id, rho_h_l2)
         write(unit_id,*)
         close(unit_id)
      end if

      previous_time = time
      previous_rho_h_integral = rho_h_integral
      previous_cumulative_advective_boundary_rho_h_flux = cumulative_advective_boundary_rho_h_flux
      previous_cumulative_conductive_boundary_flux = cumulative_conductive_boundary_flux
      previous_cumulative_qrad_integral = cumulative_qrad_integral_total
      previous_cumulative_rho_species_hdiff_integral = cumulative_rho_species_hdiff_integral_total
      previous_operator_consistent_rho_h_integral = operator_consistent_rho_h_integral
      previous_cumulative_energy_update_delta_integral = energy%cumulative_energy_update_delta_integral
      previous_cumulative_energy_update_rhs_integral = energy%cumulative_energy_update_rhs_integral
      previous_initialized = .true.

   end subroutine write_enthalpy_energy_budget_diagnostics








end module mod_output
