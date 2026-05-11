!> Basic mesh data types and finalization.
!!
!! This module defines the core geometry structures (cells, faces, patches, mesh)
!! used throughout the solver.
module mod_mesh_types
   use mod_kinds, only : rk, zero, name_len
   implicit none

   private

   public :: cell_t, face_t, patch_t, mesh_t
   public :: mesh_finalize

   integer, parameter, public :: max_cell_faces = 6

   !> Hexahedral cell data structure.
   type :: cell_t
      integer :: id = 0            !< Unique cell ID
      integer :: nodes(8) = 0      !< Indices of the 8 nodes defining the hex
      real(rk) :: center(3) = zero !< Cell centroid coordinates
      real(rk) :: volume = zero    !< Cell volume
   end type cell_t

   !> Mesh face data structure.
   type :: face_t
      integer :: id = 0                !< Unique face ID
      integer :: owner = 0             !< ID of the owner cell
      integer :: neighbor = 0          !< ID of the neighbor cell (0 if boundary)
      integer :: patch = 0             !< ID of the patch (if boundary)
      integer :: periodic_face = 0     !< ID of the corresponding periodic face
      integer :: periodic_neighbor = 0 !< ID of the neighbor across a periodic boundary
      real(rk) :: normal(3) = zero     !< Face normal vector (outward from owner)
      real(rk) :: area = zero          !< Face area
      real(rk) :: center(3) = zero     !< Face centroid coordinates
   end type face_t

   !> Boundary patch data structure.
   type :: patch_t
      integer :: id = 0                     !< Patch ID
      character(len=name_len) :: name = ""  !< Patch name (e.g., "inlet", "wall")
      integer :: nfaces = 0                 !< Number of faces in the patch
      integer, allocatable :: face_ids(:)   !< List of face IDs belonging to this patch
   end type patch_t

   !> Top-level mesh container.
   type :: mesh_t
      integer :: npoints = 0  !< Total number of points
      integer :: ncells = 0   !< Total number of cells
      integer :: nfaces = 0   !< Total number of faces
      integer :: npatches = 0 !< Total number of boundary patches

      real(rk), allocatable :: points(:,:)     !< Node coordinates (3, npoints)
      type(cell_t), allocatable :: cells(:)    !< Array of all cells
      type(face_t), allocatable :: faces(:)    !< Array of all faces
      type(patch_t), allocatable :: patches(:) !< Array of all patches

      integer, allocatable :: ncell_faces(:)   !< Number of faces per cell
      integer, allocatable :: cell_faces(:,:)  !< IDs of faces for each cell (max_cell_faces, ncells)
   end type mesh_t

contains

   !> Deallocate all mesh-related arrays.
   !!
   !! @param m Mesh data structure to finalize.
   subroutine mesh_finalize(m)
      type(mesh_t), intent(inout) :: m

      integer :: p

      if (allocated(m%points)) deallocate(m%points)
      if (allocated(m%cells)) deallocate(m%cells)
      if (allocated(m%faces)) deallocate(m%faces)
      if (allocated(m%ncell_faces)) deallocate(m%ncell_faces)
      if (allocated(m%cell_faces)) deallocate(m%cell_faces)

      if (allocated(m%patches)) then
         do p = 1, size(m%patches)
            if (allocated(m%patches(p)%face_ids)) deallocate(m%patches(p)%face_ids)
         end do
         deallocate(m%patches)
      end if

      m%npoints = 0
      m%ncells = 0
      m%nfaces = 0
      m%npatches = 0
   end subroutine mesh_finalize

end module mod_mesh_types

