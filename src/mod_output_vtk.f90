module mod_output
   use mod_kinds, only : rk, zero, path_len, fatal_error
   use mod_input, only : case_params_t
   use mod_mesh_types, only : mesh_t
   use mod_mpi_flow, only : flow_mpi_t
   use mod_fields, only : flow_fields_t
   use mod_flow_projection, only : solver_stats_t
   implicit none

   private

   public :: prepare_output, write_diagnostics_header, write_diagnostics_row
   public :: write_vtu_unstructured, write_mesh_summary, write_pvd_collection

contains

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


   subroutine write_diagnostics_header(params, flow)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow

      integer :: unit_id
      character(len=path_len + 32) :: filename

      if (flow%rank /= 0 .or. .not. params%write_diagnostics) return

      filename = trim(params%output_dir)//'/diagnostics.csv'
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a)') 'step,time,dt,max_divergence,rms_divergence,net_boundary_flux,kinetic_energy,pressure_iterations,pressure_residual'

      close(unit_id)
   end subroutine write_diagnostics_header


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

      write(unit_id,'(i0,",",es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8,",",es16.8,",",i0,",",es16.8)') &
         step, time, params%dt, stats%max_divergence, stats%rms_divergence, &
         stats%net_boundary_flux, stats%kinetic_energy, stats%pressure_iterations, &
         stats%pressure_residual

      close(unit_id)
   end subroutine write_diagnostics_row


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


   subroutine write_vtu_unstructured(params, flow, mesh, fields, step)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow
      type(mesh_t), intent(in) :: mesh
      type(flow_fields_t), intent(in) :: fields
      integer, intent(in) :: step

      integer :: unit_id
      integer :: p, c, n
      character(len=path_len + 64) :: filename

      if (flow%rank /= 0 .or. .not. params%write_vtu) return

      call validate_hex_connectivity(mesh)

      ! Keep the old subroutine name so main.f90 does not need changes.
      ! This now writes legacy ASCII VTK, not XML VTU.
      write(filename,'(a,"/flow_",i6.6,".vtk")') trim(params%output_dir), step
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a)') '# vtk DataFile Version 3.0'
      write(unit_id,'(a,i0)') 'lowmach_react_hex step ', step
      write(unit_id,'(a)') 'ASCII'
      write(unit_id,'(a)') 'DATASET UNSTRUCTURED_GRID'

      ! ------------------------------------------------------------
      ! Points
      ! ------------------------------------------------------------
      write(unit_id,'(a,i0,a)') 'POINTS ', mesh%npoints, ' double'

      do p = 1, mesh%npoints
         write(unit_id,'(3(es24.16,1x))') mesh%points(1,p), mesh%points(2,p), mesh%points(3,p)
      end do

      ! ------------------------------------------------------------
      ! Cells
      !
      ! Legacy VTK format:
      ! CELLS number_of_cells total_integer_count
      !
      ! Each hexahedron line has:
      ! 8 node0 node1 node2 node3 node4 node5 node6 node7
      !
      ! So each cell contributes 9 integers.
      ! ------------------------------------------------------------
      write(unit_id,'(a,i0,1x,i0)') 'CELLS ', mesh%ncells, 9 * mesh%ncells

      do c = 1, mesh%ncells
         write(unit_id,'(i0,1x,8(i0,1x))') 8, (mesh%cells(c)%nodes(n) - 1, n = 1, 8)
      end do

      ! ------------------------------------------------------------
      ! Cell types
      !
      ! VTK_HEXAHEDRON = 12
      ! ------------------------------------------------------------
      write(unit_id,'(a,i0)') 'CELL_TYPES ', mesh%ncells

      do c = 1, mesh%ncells
         write(unit_id,'(i0)') 12
      end do

      ! ------------------------------------------------------------
      ! Cell data
      ! ------------------------------------------------------------
      write(unit_id,'(a,i0)') 'CELL_DATA ', mesh%ncells

      write(unit_id,'(a)') 'VECTORS velocity double'
      do c = 1, mesh%ncells
         write(unit_id,'(3(es24.16,1x))') fields%u(1,c), fields%u(2,c), fields%u(3,c)
      end do

      call write_vtk_cell_scalar(unit_id, 'pressure', fields%p)
      call write_vtk_cell_scalar(unit_id, 'divergence', fields%div)

      close(unit_id)
   end subroutine write_vtu_unstructured


   subroutine write_pvd_collection(params, flow, nsteps, output_interval, dt)
      type(case_params_t), intent(in) :: params
      type(flow_mpi_t), intent(in) :: flow
      integer, intent(in) :: nsteps
      integer, intent(in) :: output_interval
      real(rk), intent(in) :: dt

      integer :: unit_id
      integer :: step
      character(len=path_len + 32) :: filename
      character(len=64) :: vtk_name

      if (flow%rank /= 0 .or. .not. params%write_vtu) return

      filename = trim(params%output_dir)//'/flow.pvd'
      open(newunit=unit_id, file=trim(filename), status='replace', action='write')

      write(unit_id,'(a)') '<?xml version="1.0"?>'
      write(unit_id,'(a)') '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian">'
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

      subroutine write_dataset_line(unit_id, step_id, time_value)
         integer, intent(in) :: unit_id
         integer, intent(in) :: step_id
         real(rk), intent(in) :: time_value

         character(len=32) :: time_text

         write(vtk_name,'("flow_",i6.6,".vtk")') step_id
         write(time_text,'(es16.8)') time_value
         time_text = adjustl(time_text)

         write(unit_id,'(a,a,a,a,a)') '    <DataSet timestep="', trim(time_text), &
            '" group="" part="0" file="', trim(vtk_name), '"/>'
      end subroutine write_dataset_line

   end subroutine write_pvd_collection


   subroutine write_vtk_cell_scalar(unit_id, name, field)
      integer, intent(in) :: unit_id
      character(len=*), intent(in) :: name
      real(rk), intent(in) :: field(:)

      integer :: i

      write(unit_id,'(a,a,a)') 'SCALARS ', trim(name), ' double 1'
      write(unit_id,'(a)') 'LOOKUP_TABLE default'

      do i = 1, size(field)
         write(unit_id,'(es24.16)') field(i)
      end do
   end subroutine write_vtk_cell_scalar


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