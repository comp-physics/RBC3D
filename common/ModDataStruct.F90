module ModDataStruct

  use ModDataTypes

  implicit none
#include "../petsc_include.h"

  private

  public :: t_Spline, &
    t_Rbc, &
    t_RbcPolarPatch, &
    t_Wall, &
    t_SourceList, &
    t_TargetList, &
    t_NbrRbcList

!**********************************************************************
! 2D periodic spline surface of periodicity 2*PI
!
! Note:
!  Can contain a linearly increasing part which has slopes kx and ky
!
!======================================================================
! m, n -- number of mesh points
! nvar -- number of variables
! hx, hy -- mesh spacing
! u -- values of variables on the mesh
! u1, u2 -- derivatives along the two directions
! u12 -- mixed derivative
! kx, ky -- slope of the linear part of u such that
!   (u - kx*x) is periodic in x-direction
!   (u - ky*y) is periodic in y-direction
!**********************************************************************
  type t_spline
    integer :: m, n, nvar
    real(WP) :: hx, hy, ihx, ihy
    real(WP),pointer :: u(:,:,:)
    real(WP),pointer :: u1(:,:,:), u2(:,:,:), u12(:,:,:)
    real(WP),pointer :: kx(:), ky(:)
  end type t_spline


!**********************************************************************
! Red blood cell surface
!
!======================================================================
! Mesh definition
!
! nlat0, nlon0 -- number of latitudinal modes and longitudinal modes
!
! nlat, nlon -- number of mesh points
!
! th -- collatitude angles
! phi -- longitudinal angles
! w -- integration weights (uniform along longitudinal direction)
!
! x(ilat,ilon,:) -- mesh point coordinate
!======================================================================
! Geometry
!
! detJ(ilat,ilon) -- Jacobian
! a3(ilat,ilon,:) -- surface normal
!
! a1(ilat,ilon,:), a2(ilat,ilon,:) -- surface tangents
! a(ilat,ilon,:,:) -- the first fundamental form
! a1_rcp(ilat,ilon,:), a2_rcp(ilat,ilon,:) -- reciprocal tangents
! a_rcp(ilat,ilon,:,:) -- inverse of a(ilat,ilon,:,:)
!
! b(ilat,ilon,:,:) -- the second fundamental form
!
! vol -- cell volume
! area -- surface area
! xc -- center of the cell
! meshSize -- representative mesh size, sqrt( area/(nlat*nlat) )
!
!======================================================================
! Surface source density
!
! f(ilat,ilon,:) -- surface force density
! g(ilat,ilon,:) -- double-layer potential density
! 
!======================================================================
! Spline interpolation
!
!  spln_x -- x
!  spln_a3 -- a3
!  spln_J -- detJ
!  spln_fJ -- f*detJ
!  spln_gJ -- g*detJ
!
!======================================================================
! Cell type
!  celltype --   1  red cell ; 2 leukocyte ; 3 platelet
!
!======================================================================
! Solid body modes
!  qq g(ilat,ilon,:,6) -- 6-vector: 3 translational, 3 rotational
!
!======================================================================
!  patch -- polar patch for singular integration
!======================================================================
! Material properties
!
! ES, ED, EB -- shear, dilatation and bending modulus
!  
!======================================================================
  type t_RBC
    integer :: nlat0, nlon0
    integer :: nlat, nlon
    real(WP),dimension(:),pointer :: th, phi, w
    real(WP),dimension(:),pointer :: wshags, wshsgs, wshses, wvhsgs
    real(WP),dimension(:,:,:),pointer :: x
    real(WP),dimension(:,:,:),pointer :: a1, a2, a1_rcp, a2_rcp, a3
    real(WP),dimension(:,:),pointer :: detj
    real(WP),dimension(:,:,:,:),pointer :: a, a_rcp, b

    real(WP),dimension(:,:,:),pointer :: f, g, v  ! single, double, velocity
    ! vel and double are the same for finite viscosity ratio

    real(WP),dimension(:,:,:,:),pointer :: qq

    real(WP) :: xc(3), vol, area, startingArea !is area here the local area value?
    real(WP) :: meshSize

    type(t_spline) :: spln_x, spln_a3, spln_detJ, spln_FdetJ, spln_GdetJ
    type(t_RbcPolarPatch),pointer :: patch

    real(WP) :: ES, ED, EB

    integer :: ID   ! the unique ID

    integer :: celltype

  end type t_RBC

!**********************************************************************
! Polar coordinate patch for a spherical surface 
!
! nlat, nlon
! radius -- patch radius
! nrad, nazm -- number of patch points along radial and 
!       azimuthal directions
! thL, phiL -- local coordiantes of patch mesh points
! w -- weight of integration weight (including masking fucntion)
! thG(:,:,i,j), phiG(:,:,i,j) -- global latitudinal and longitudinal
!       coordinates of a patch centered at (i,j)-th mesh point
  type t_RbcPolarPatch
    integer :: nlat, nlon
    real(WP) :: radius

    integer :: nrad, nazm
    real(WP),dimension(:),pointer :: thL, phiL, w;
    real(WP),dimension(:,:,:,:),pointer :: thG, phiG;
  end type t_RbcPolarPatch

!**********************************************************************
! Surface representation by a boundary-element format
!
! Members:
! nvert, nele -- number of vertices and elements
!
! x(i,:) -- current and equilibrium positions of points
! f -- force density at each vertex
! g -- velocity
!
! e2v(i,:) -- element to vertex
! v2v(i) -- for periodic boundary conditions, pointing to the same vertex
!           which is shifted by one period
! Constraint: v2v(i) <= i v2v(v2v(i)) = v2v(i)
!
! indxVertGlb -- global vertex index
!
! a3 -- element normal
! area -- element area
! epsDist -- threshhold distance for singular integration
! areaTot -- total surface area
!
! lhs -- lhs matrix for no-slip boundary condition
! 
!**********************************************************************
  type t_Wall
    integer :: nvert, nele
    real(WP),pointer :: x(:,:), f(:,:), g(:,:)
    integer,pointer :: e2v(:,:), v2v(:)

    integer,pointer :: indxVertGlb(:)

    real(WP),pointer :: a3(:,:), area(:), epsDist(:)
    real(WP) :: areaTot

    Mat :: lhs

    integer :: ID
  end type t_Wall


!**********************************************************************
! Collection of singularity sources
!  nPoint  -- number of points
!  x(i,:) -- coordinate of the i-th point
!  f(i,:) -- force singularity at the i-th point
!  g(i,:) -- double-layer potential source singularity at the i-th piont
!  a3(i,:) -- normal direction
!  lam(i) -- viscRatio
!  
!  indx(:,0:2) -- index
!     for RBC surface, (surfId, ilat, ilon)
!     for wall surface, (surfId, iele, -1)
!  
!  Nc -- number of cells
!  iLbNc -- Nc/Lb, for indentifying which cell a point lies in
!  hoc -- first point in the cell
!  next -- next point
  type t_SourceList
    integer :: nPoint
    real(WP),pointer,dimension(:,:) :: x, f, g, a3
    real(WP),pointer,dimension(:) :: lam
    real(WP),pointer,dimension(:) :: Bcoef
    integer,pointer :: indx(:,:)

    integer :: Nc(3)
    real(WP) :: iLbNc(3)
    integer,pointer :: hoc(:,:,:), next(:)
  end type t_SourceList

!**********************************************************************
! Collection of target points
!  nPoint -- number of points
!  x -- coordinate and velocities of points
!  indx -- same as those in t_SourceList
!  active -- wheather the target point is active
  type t_TargetList
    integer :: nPoint
    real(WP),pointer,dimension(:,:) :: x
    real(WP),pointer,dimension(:) :: lam
    real(WP),pointer,dimension(:) :: Acoef
    integer,pointer :: indx(:,:)
    logical,pointer :: active(:)
  end type t_TargetList

!**********************************************************************
! Collection of a cells that are too close to a point
!  N -- number of cells in the list
!  dist -- smallest distance to the cell
!  indx -- index of the mesh point whose distance to the target point gives "dist"
!    indx(0) -- surfId
!    indx(1:2) -- mesh index on the surface
  type t_NbrRbcList
    integer :: N
    real(WP),pointer :: dist(:)
    integer,pointer :: indx(:,:)
  end type t_NbrRbcList

!**********************************************************************

end module ModDataStruct
