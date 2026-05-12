!> Mesh data structures and geometry definitions for hexahedral grids.
!!
!! This module defines the core topological and geometric entities used by the 
!! finite-volume solver. The solver is specifically designed for 
!! **hexahedral cells** (8 nodes, 6 faces), and these structures encapsulate 
!! all necessary connectivity and metric information (volumes, areas, normals).
!!
!! The mesh hierarchy is structured to support efficient flux-based 
!! computations:
!!
!! 1. **`cell_t`**: Holds volumetric data and node connectivity. It represents 
!!    the primary control volume.
!! 2. **`face_t`**: The fundamental entity for flux evaluation. Manages 
!!    owner/neighbor relationships and boundary associations.
!! 3. **`patch_t`**: Groups faces into named boundary regions (e.g., "Inlet", 
!!    "Outlet") for application of Boundary Conditions.
!! 4. **`mesh_t`**: The global container for all nodes, cells, and faces.
!!
!! **Gmsh Hexahedron Node Ordering (8-node):**
!! ```
!!       7 --------- 6
!!      /|          /|
!!     4 --------- 5 |
!!     | |         | |
!!     | 3 --------| 2
!!     |/          |/
!!     0 --------- 1
!! ```
module mod_mesh_types
   use mod_kinds, only : rk, zero, name_len
   implicit none

   private

   public :: cell_t, face_t, patch_t, mesh_t
   public :: mesh_finalize

   !> Maximum number of faces per cell. 
   !! Fixed at 6 for hexahedral meshes.
   integer, parameter, public :: max_cell_faces = 6

   !> Hexahedral cell data structure.
   !! Contains local geometric properties and node-based connectivity.
   type :: cell_t
      integer :: id = 0            !< Unique global cell index.
      integer :: nodes(8) = 0      !< Indices of the 8 nodes defining the hex (standard Gmsh ordering).
      real(rk) :: center(3) = zero !< Coordinates (x, y, z) of the cell centroid.
      real(rk) :: volume = zero    !< Total cell volume [m^3].
   end type cell_t

   !> Mesh face data structure.
   !! Faces are the primary entities for flux calculation. Each face is 
   !! associated with an 'owner' and potentially a 'neighbor' cell.
   type :: face_t
      integer :: id = 0                !< Unique global face index.
      integer :: owner = 0             !< ID of the owner cell (always exists).
      integer :: neighbor = 0          !< ID of the neighbor cell (0 for boundary faces).
      integer :: patch = 0             !< ID of the associated boundary patch (0 for internal faces).
      integer :: periodic_face = 0     !< ID of the matching face on the opposite periodic boundary.
      integer :: periodic_neighbor = 0 !< ID of the neighbor cell accessed through periodic matching.
      real(rk) :: normal(3) = zero     !< Unit normal vector pointing outward from the owner cell.
      real(rk) :: area = zero          !< Surface area of the face [m^2].
      real(rk) :: center(3) = zero     !< Coordinates (x, y, z) of the face centroid.
   end type face_t

   !> Boundary patch data structure.
   !! Patches group faces together to allow bulk application of boundary 
   !! conditions (e.g., all faces in the "inlet" patch).
   type :: patch_t
      integer :: id = 0                     !< Unique patch index.
      character(len=name_len) :: name = ""  !< Descriptive name (e.g., "inlet", "wall", "outlet").
      integer :: nfaces = 0                 !< Total number of faces assigned to this patch.
      integer, allocatable :: face_ids(:)   !< List of global face IDs belonging to this patch.
   end type patch_t

   !> Top-level mesh container.
   !! This structure holds all geometric data and connectivity arrays for 
   !! the entire computational domain.
   type :: mesh_t
      integer :: ncells = 0   !< Total number of hexahedral cells.
      integer :: nfaces = 0   !< Total number of faces (internal + boundary).
      integer :: npoints = 0  !< Total number of unique nodes in the mesh.
      integer :: npatches = 0 !< Total number of boundary regions defined.

      real(rk), allocatable :: points(:,:)     !< Node coordinates (3, npoints).
      type(cell_t), allocatable :: cells(:)    !< Array containing all cell geometric data.
      type(face_t), allocatable :: faces(:)    !< Array containing all face connectivity data.
      type(patch_t), allocatable :: patches(:) !< Array of all boundary patches.

      integer, allocatable :: ncell_faces(:)   !< Number of faces belonging to each cell (usually 6).
      integer, allocatable :: cell_faces(:,:)  !< Indices of the faces for each cell (max_cell_faces, ncells).
   end type mesh_t

contains

   !> Safely deallocates all heap-allocated arrays in the mesh structure.
   !!
   !! This routine should be called during solver shutdown or mesh reload 
   !! to prevent memory leaks. It nullifies all counts after deallocation.
   !!
   !! @param m The mesh object to be cleared.
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

