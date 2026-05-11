module mod_mesh_io
   use mod_kinds, only : rk, path_len, fatal_error
   use mod_mesh_types, only : mesh_t, mesh_finalize, max_cell_faces
   implicit none

   private

   public :: read_native_mesh

contains

   subroutine read_native_mesh(mesh_dir, m)
      character(len=*), intent(in) :: mesh_dir
      type(mesh_t), intent(inout) :: m

      call mesh_finalize(m)
      call read_points(trim(mesh_dir)//'/points.dat', m)
      call read_cells(trim(mesh_dir)//'/cells.dat', m)
      call read_faces(trim(mesh_dir)//'/faces.dat', m)
      call read_patches(trim(mesh_dir)//'/patches.dat', m)
      call read_periodic_optional(trim(mesh_dir)//'/periodic.dat', m)
      call build_cell_faces(m)
   end subroutine read_native_mesh


   subroutine read_points(filename, m)
      character(len=*), intent(in) :: filename
      type(mesh_t), intent(inout) :: m

      integer :: unit_id, ios, i, id
      real(rk) :: x, y, z

      open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) call fatal_error('mesh_io', 'could not open '//trim(filename))

      read(unit_id, *, iostat=ios) m%npoints
      if (ios /= 0 .or. m%npoints <= 0) call fatal_error('mesh_io', 'invalid points header')

      allocate(m%points(3, m%npoints))
      do i = 1, m%npoints
         read(unit_id, *, iostat=ios) id, x, y, z
         if (ios /= 0) call fatal_error('mesh_io', 'failed reading point')
         if (id < 1 .or. id > m%npoints) call fatal_error('mesh_io', 'point id out of range')
         m%points(:, id) = [x, y, z]
      end do

      close(unit_id)
   end subroutine read_points


   subroutine read_cells(filename, m)
      character(len=*), intent(in) :: filename
      type(mesh_t), intent(inout) :: m

      integer :: unit_id, ios, i, id
      integer :: nodes(8)
      real(rk) :: cx, cy, cz, volume

      open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) call fatal_error('mesh_io', 'could not open '//trim(filename))

      read(unit_id, *, iostat=ios) m%ncells
      if (ios /= 0 .or. m%ncells <= 0) call fatal_error('mesh_io', 'invalid cells header')

      allocate(m%cells(m%ncells))
      do i = 1, m%ncells
         read(unit_id, *, iostat=ios) id, nodes, cx, cy, cz, volume
         if (ios /= 0) call fatal_error('mesh_io', 'failed reading cell')
         if (id < 1 .or. id > m%ncells) call fatal_error('mesh_io', 'cell id out of range')
         if (volume <= 0.0_rk) call fatal_error('mesh_io', 'cell volume must be positive')
         m%cells(id)%id = id
         m%cells(id)%nodes = nodes
         m%cells(id)%center = [cx, cy, cz]
         m%cells(id)%volume = volume
      end do

      close(unit_id)
   end subroutine read_cells


   subroutine read_faces(filename, m)
      character(len=*), intent(in) :: filename
      type(mesh_t), intent(inout) :: m

      integer :: unit_id, ios, i, id
      integer :: owner, neighbor, patch
      real(rk) :: nx, ny, nz, area, cx, cy, cz

      open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) call fatal_error('mesh_io', 'could not open '//trim(filename))

      read(unit_id, *, iostat=ios) m%nfaces
      if (ios /= 0 .or. m%nfaces <= 0) call fatal_error('mesh_io', 'invalid faces header')

      allocate(m%faces(m%nfaces))
      do i = 1, m%nfaces
         read(unit_id, *, iostat=ios) id, owner, neighbor, patch, nx, ny, nz, area, cx, cy, cz
         if (ios /= 0) call fatal_error('mesh_io', 'failed reading face')
         if (id < 1 .or. id > m%nfaces) call fatal_error('mesh_io', 'face id out of range')
         if (owner < 1 .or. owner > m%ncells) call fatal_error('mesh_io', 'face owner out of range')
         if (neighbor < 0 .or. neighbor > m%ncells) call fatal_error('mesh_io', 'face neighbor out of range')
         if (area <= 0.0_rk) call fatal_error('mesh_io', 'face area must be positive')
         m%faces(id)%id = id
         m%faces(id)%owner = owner
         m%faces(id)%neighbor = neighbor
         m%faces(id)%patch = patch
         m%faces(id)%normal = [nx, ny, nz]
         m%faces(id)%area = area
         m%faces(id)%center = [cx, cy, cz]
      end do

      close(unit_id)
   end subroutine read_faces


   subroutine read_patches(filename, m)
      character(len=*), intent(in) :: filename
      type(mesh_t), intent(inout) :: m

      integer :: unit_id, ios, p, id, nfaces
      character(len=64) :: name

      open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) call fatal_error('mesh_io', 'could not open '//trim(filename))

      read(unit_id, *, iostat=ios) m%npatches
      if (ios /= 0 .or. m%npatches < 0) call fatal_error('mesh_io', 'invalid patches header')

      allocate(m%patches(m%npatches))
      do p = 1, m%npatches
         read(unit_id, *, iostat=ios) id, name, nfaces
         if (ios /= 0) call fatal_error('mesh_io', 'failed reading patch header')
         if (id < 1 .or. id > m%npatches) call fatal_error('mesh_io', 'patch id out of range')
         m%patches(id)%id = id
         m%patches(id)%name = name
         m%patches(id)%nfaces = nfaces
         allocate(m%patches(id)%face_ids(nfaces))
         if (nfaces > 0) then
            read(unit_id, *, iostat=ios) m%patches(id)%face_ids
            if (ios /= 0) call fatal_error('mesh_io', 'failed reading patch face ids')
         else
            read(unit_id, *, iostat=ios)
         end if
      end do

      close(unit_id)
   end subroutine read_patches


   subroutine read_periodic_optional(filename, m)
      character(len=*), intent(in) :: filename
      type(mesh_t), intent(inout) :: m

      integer :: unit_id, ios, nlinks, i
      integer :: face_id, pair_face_id, neighbor_cell

      open(newunit=unit_id, file=trim(filename), status='old', action='read', iostat=ios)
      if (ios /= 0) return

      read(unit_id, *, iostat=ios) nlinks
      if (ios /= 0 .or. nlinks < 0) call fatal_error('mesh_io', 'invalid periodic header')

      do i = 1, nlinks
         read(unit_id, *, iostat=ios) face_id, pair_face_id, neighbor_cell
         if (ios /= 0) call fatal_error('mesh_io', 'failed reading periodic link')
         if (face_id < 1 .or. face_id > m%nfaces) call fatal_error('mesh_io', 'periodic face id out of range')
         if (pair_face_id < 1 .or. pair_face_id > m%nfaces) call fatal_error('mesh_io', 'periodic pair face id out of range')
         if (neighbor_cell < 1 .or. neighbor_cell > m%ncells) call fatal_error('mesh_io', 'periodic neighbor cell out of range')
         m%faces(face_id)%periodic_face = pair_face_id
         m%faces(face_id)%periodic_neighbor = neighbor_cell
      end do

      close(unit_id)
   end subroutine read_periodic_optional


   subroutine build_cell_faces(m)
      type(mesh_t), intent(inout) :: m

      integer :: f, c

      allocate(m%ncell_faces(m%ncells))
      allocate(m%cell_faces(max_cell_faces, m%ncells))
      m%ncell_faces = 0
      m%cell_faces = 0

      do f = 1, m%nfaces
         c = m%faces(f)%owner
         call append_cell_face(m, c, f)

         c = m%faces(f)%neighbor
         if (c > 0) call append_cell_face(m, c, f)
      end do

      do c = 1, m%ncells
         if (m%ncell_faces(c) /= max_cell_faces) then
            call fatal_error('mesh_io', 'each v1 cuboid cell must have exactly six faces')
         end if
      end do
   end subroutine build_cell_faces


   subroutine append_cell_face(m, cell_id, face_id)
      type(mesh_t), intent(inout) :: m
      integer, intent(in) :: cell_id
      integer, intent(in) :: face_id

      if (m%ncell_faces(cell_id) >= max_cell_faces) then
         call fatal_error('mesh_io', 'cell has more than six faces')
      end if

      m%ncell_faces(cell_id) = m%ncell_faces(cell_id) + 1
      m%cell_faces(m%ncell_faces(cell_id), cell_id) = face_id
   end subroutine append_cell_face

end module mod_mesh_io

