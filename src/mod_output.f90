!> Output management for VTK visualization and diagnostics (XML VTU format).
!!
!! This module handles the generation of simulation results in modern XML-based 
!! VTK format (.vtu) and CSV-based global diagnostics. It manages the 
!! creation of the output directory, writing mesh summaries, and 
!! generating PVD collection files for time-series visualization in ParaView.
module mod_output
   use mod_kinds, only : rk, zero, path_len, fatal_error
   use mod_input, only : case_params_t
   use mod_mesh_types, only : mesh_t
   use mod_mpi_flow, only : flow_mpi_t, flow_gather_owned_scalar_root, &
                            flow_gather_owned_matrix_root
   use mod_fields, only : flow_fields_t
   use mod_flow_projection, only : solver_stats_t
   use mod_species, only : species_fields_t
   implicit none

   private

   public :: prepare_output, write_diagnostics_header, write_diagnostics_row
   public :: write_vtu_unstructured, write_mesh_summary, write_pvd_collection

contains

   !> Creates the output directory specified in the case parameters.
   subroutine prepare_output(params, flow)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow

      integer :: exitstat
      character(len=path_len + 16) :: command

      if (flow%rank /= 0) return

      command = 'mkdir -p '//trim(params%output_dir)
      call execute_command_line(trim(command), exitstat=exitstat)

      if (exitstat /= 0) then
         call fatal_error('output', 'failed to create output directory')
      end if
   end subroutine prepare_output


   !> Writes the CSV header for global simulation diagnostics.
   subroutine write_diagnostics_header(params, flow)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow

      integer :: unit_id
      character(len=path_len + 32) :: filename

      if (flow%rank /= 0 .or. .not. params%write_diagnostics) return

      filename = trim(params%output_dir)//'/diagnostics.csv'
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

      filename = trim(params%output_dir)//'/diagnostics.csv'
      open(newunit=unit_id, file=trim(filename), status='old', position='append', action='write')

      write(unit_id,'(i0,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,'// &
                    'i0,a,i0,a,i0,a,es16.8,a,i0,a,es16.8,a,es16.8,a,es16.8,a,'// &
                    'es16.8,a,es16.8,a,es16.8)') &
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
   !! - Divergence (Cell Scalar)
   !! - Species Mass Fractions (Cell Scalars, if enabled)
   subroutine write_vtu_unstructured(params, flow, mesh, fields, species, step)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(inout) :: flow
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(in) :: fields
      type(species_fields_t), intent(in) :: species
      integer, intent(in) :: step

      integer :: unit_id
      integer :: p, c, n, k
      character(len=path_len + 64) :: filename
      real(rk), allocatable :: out_u(:,:), out_p(:), out_div(:), out_Y(:,:)

      if (.not. params%write_vtu) return

      allocate(out_u(3, mesh%ncells))
      allocate(out_p(mesh%ncells))
      allocate(out_div(mesh%ncells))

      call flow_gather_owned_matrix_root(flow, fields%u, out_u)
      call flow_gather_owned_scalar_root(flow, fields%p, out_p)
      call flow_gather_owned_scalar_root(flow, fields%div, out_div)

      if (params%enable_species .and. params%nspecies > 0) then
         allocate(out_Y(params%nspecies, mesh%ncells))
         call flow_gather_owned_matrix_root(flow, species%Y, out_Y)
      end if

      if (flow%rank /= 0) then
         if (allocated(out_u)) deallocate(out_u)
         if (allocated(out_p)) deallocate(out_p)
         if (allocated(out_div)) deallocate(out_div)
         if (allocated(out_Y)) deallocate(out_Y)
         return
      end if

      call validate_hex_connectivity(mesh)

      write(filename,'(a,"/flow_",i6.6,".vtu")') trim(params%output_dir), step
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a)') '<?xml version="1.0"?>'
      write(unit_id,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
      write(unit_id,'(a)') '  <UnstructuredGrid>'
      write(unit_id,'(a,i0,a,i0,a)') '    <Piece NumberOfPoints="', mesh%npoints, &
                                      '" NumberOfCells="', mesh%ncells, '">'

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
      do c = 1, mesh%ncells
         write(unit_id,'(3(es24.16,1x))') out_u(1,c), out_u(2,c), out_u(3,c)
      end do
      write(unit_id,'(a)') '        </DataArray>'

      call write_vtu_cell_scalar(unit_id, 'pressure', out_p)
      call write_vtu_cell_scalar(unit_id, 'divergence', out_div)

      if (params%enable_species .and. params%nspecies > 0) then
         do k = 1, params%nspecies
            call write_vtu_cell_scalar(unit_id, 'Y_'//trim(params%species_name(k)), out_Y(k,:))
         end do

         write(unit_id,'(a)') '        <DataArray type="Float64" Name="sum_Y" format="ascii">'
         do c = 1, mesh%ncells
            write(unit_id,'(es24.16)') sum(out_Y(:,c))
         end do
         write(unit_id,'(a)') '        </DataArray>'
      end if

      ! Helpful debug scalar so you can color by cell_id in ParaView.
      write(unit_id,'(a)') '        <DataArray type="Int32" Name="cell_id" format="ascii">'
      do c = 1, mesh%ncells
         write(unit_id,'(i0)') c
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a)') '      </CellData>'

      write(unit_id,'(a)') '      <Points>'
      write(unit_id,'(a)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
      do p = 1, mesh%npoints
         write(unit_id,'(3(es24.16,1x))') mesh%points(1,p), mesh%points(2,p), mesh%points(3,p)
      end do
      write(unit_id,'(a)') '        </DataArray>'
      write(unit_id,'(a)') '      </Points>'

      write(unit_id,'(a)') '      <Cells>'

      write(unit_id,'(a)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
      do c = 1, mesh%ncells
         write(unit_id,'(8(i0,1x))') (mesh%cells(c)%nodes(n) - 1, n = 1, 8)
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
      do c = 1, mesh%ncells
         write(unit_id,'(i0)') 8 * c
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a)') '        <DataArray type="UInt8" Name="types" format="ascii">'
      do c = 1, mesh%ncells
         write(unit_id,'(i0)') 12
      end do
      write(unit_id,'(a)') '        </DataArray>'

      write(unit_id,'(a)') '      </Cells>'

      write(unit_id,'(a)') '    </Piece>'
      write(unit_id,'(a)') '  </UnstructuredGrid>'
      write(unit_id,'(a)') '</VTKFile>'

      close(unit_id)

      deallocate(out_u)
      deallocate(out_p)
      deallocate(out_div)
      if (allocated(out_Y)) deallocate(out_Y)
   end subroutine write_vtu_unstructured


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

      filename = trim(params%output_dir)//'/flow.pvd'
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

         write(vtu_name,'("flow_",i6.6,".vtu")') step_id
         write(time_text,'(es16.8)') time_value
         time_text = adjustl(time_text)

         write(unit_id,'(a,a,a,a,a)') '    <DataSet timestep="', trim(time_text), &
            '" group="" part="0" file="', trim(vtu_name), '"/>'
      end subroutine write_dataset_line

   end subroutine write_pvd_collection


   !> Internal helper to write a scalar field to a VTU file.
   subroutine write_vtu_cell_scalar(unit_id, name, field)
      integer, intent(in) :: unit_id
      character(len=*), intent(in) :: name
      real(rk), intent(in) :: field(:)

      integer :: i

      write(unit_id,'(a,a,a)') '        <DataArray type="Float64" Name="', trim(name), '" format="ascii">'

      do i = 1, size(field)
         write(unit_id,'(es24.16)') field(i)
      end do

      write(unit_id,'(a)') '        </DataArray>'
   end subroutine write_vtu_cell_scalar


   !> Performs sanity checks on hex connectivity before writing output.
   subroutine validate_hex_connectivity(mesh)
      type(mesh_t), intent(in) :: mesh

      integer :: c, n, node_id

      if (mesh%npoints <= 0) then
         call fatal_error('output', 'mesh has no points')
      end if

      if (mesh%ncells <= 0) then
         call fatal_error('output', 'mesh has no cells')
      end if

      do c = 1, mesh%ncells
         do n = 1, 8
            node_id = mesh%cells(c)%nodes(n)

            if (node_id < 1 .or. node_id > mesh%npoints) then
               call fatal_error('output', 'cell connectivity has node id outside 1..npoints')
            end if
         end do
      end do
   end subroutine validate_hex_connectivity

end module mod_output
