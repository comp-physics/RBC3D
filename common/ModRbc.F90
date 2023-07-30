module ModRbc

  use ModDataTypes
  use ModBasicMath
  use ModDataStruct
  use ModSphpk
  use ModSpline
  use ModConf

  implicit none

  private

  public :: RBC_Create, &
    RBC_Destroy, &
    RBC_MakeSphere, &
    RBC_MakeEllipsoid, &
    RBC_MakeLeukocyte, &
    RBC_MakePlatelet, &
    RBC_MakeBiConcave, &
    RBC_ComputeGeometry, &
    RBC_CovarGrad_Vec, &
    RBC_CovarGrad_Tensor, &
    RBC_SphProject, &
    RBC_Integral, &
    RBC_BuildSurfaceSource, &
    RBC_MakeUncurvedSickleSpheroid

  public :: Shell_ResForce
  private :: Shell_ElasStrs, &
    Shell_Bend


contains

!**********************************************************************
! Create a red blood cell surface mesh
! Arguments:
!  cell -- 
!  nlat0 -- number of latitudinal modes
!  dealias_fac -- dealiasing factor
  subroutine RBC_Create(cell, nlat0, dealias_fac)
    type(t_RBC) :: cell
    integer :: nlat0
    integer, optional :: dealias_fac

    integer :: nlat, nlon, i, j
    ! Working variables for spherepack subroutine gaqd(...)
    double precision :: work
    integer :: lwork, ierr

    ! Set up resolution
    cell%nlat0 = nlat0
    cell%nlon0 = 2*cell%nlat0
    if (present(dealias_fac)) then
       if (dealias_fac.eq.100) then
          
          cell%nlat = nlat0 + 2
          cell%nlon = 2*cell%nlat 

      else if (dealias_fac.gt.1) then

          cell%nlat = dealias_fac*nlat0
          cell%nlon = 2*cell%nlat 

       else if (dealias_fac.eq.1) then
          
          cell%nlat = nlat0 
          cell%nlon = 2*cell%nlat 
          
       end if
    else 
       cell%nlat = 3*nlat0
       cell%nlon = 2*cell%nlat
    end if

    ! Geometry mesh
    nlat = cell%nlat
    nlon = cell%nlon

    allocate(cell%th(nlat), cell%phi(nlon), cell%w(nlat) )
    allocate(cell%x(nlat,nlon,3) )
    allocate(cell%a1(nlat,nlon,3), cell%a2(nlat,nlon,3), &
            cell%a1_rcp(nlat,nlon,3), cell%a2_rcp(nlat,nlon,3), &
        cell%a3(nlat,nlon,3) )
    allocate(cell%detj(nlat,nlon) )
    allocate(cell%a(nlat,nlon,2,2), cell%a_rcp(nlat,nlon,2,2), &
            cell%b(nlat,nlon,2,2) )
    allocate(cell%f(nlat,nlon,3), cell%g(nlat,nlon,3), cell%v(nlat,nlon,3)  )

    allocate(cell%qq(nlat,nlon,3,6))

    call gaqd(nlat, cell%th, cell%w, work, lwork, ierr)
    cell%w = (TWO_PI/nlon)*cell%w
    cell%phi = (/ ((j-1)*TWO_PI/nlon, j=1,nlon) /)

    call ShAnalGau_Init(nlat, nlon, cell%wshags)
    call ShSynthGau_Init(nlat, nlon, cell%wshsgs)
    call ShSynthEqu_Init(nlat+1, nlon, cell%wshses)
    call ShGradGau_Init(nlat, nlon, cell%wvhsgs)

    call Spline_Create(cell%spln_x, nlat*2, nlon, 3)
    call Spline_Create(cell%spln_a3, nlat*2, nlon, 3)
    call spline_Create(cell%spln_detJ, nlat*2, nlon, 1)
    call Spline_Create(cell%spln_FdetJ, nlat*2, nlon, 3)
    call Spline_Create(cell%spln_GdetJ, nlat*2, nlon, 3)

  end subroutine RBC_Create

!**********************************************************************
! Destroy a red blood cell
  subroutine RBC_Destroy(cell)
    type(t_RBC) :: cell

    deallocate(cell%th, cell%phi, cell%w, cell%x)
    deallocate(cell%a1, cell%a2, cell%a1_rcp, cell%a2_rcp, cell%a3)
    deallocate(cell%detj, cell%a, cell%a_rcp, cell%b, cell%f, cell%g, cell%v, cell%qq)
    deallocate(cell%wshags, cell%wshsgs, cell%wshses, cell%wvhsgs)

    call Spline_Destroy(cell%spln_x)
    call Spline_Destroy(cell%spln_a3)
    call Spline_Destroy(cell%spln_detJ)
    call Spline_Destroy(cell%spln_FdetJ)
    call Spline_Destroy(cell%spln_GdetJ)


  end subroutine RBC_Destroy

!**********************************************************************
! Make a spherical shaped cell for modeling a larger leukocyte --- r factor is hard coded here (for now at least)
! Arguments:
!  cell -- red blood cell
!  r -- spherical radius
!  xc -- center of the sphere

  subroutine RBC_MakeLeukocyte(cell, r, xc)
    type(t_RBC) :: cell
    real(WP) :: r
    real(WP),optional :: xc(3)
  
    real(WP) :: rleuk

    rleuk = 1.4*r   ! Leukocyte expansion factor
    call RBC_MakeSphere(cell, rleuk, xc)

  end subroutine RBC_MakeLeukocyte

!**********************************************************************
! Make a spherical shaped cell
! Arguments:
!  cell -- red blood cell
!  r -- spherical radius
!  xc -- center of the sphere
  subroutine RBC_MakeSphere(cell, r, xc)
    type(t_RBC) :: cell
    real(WP) :: r
    real(WP),optional :: xc(3)

    integer :: ilat, ilon, ii
    real(WP) :: th, phi

    do ilon = 1, cell%nlon
    do ilat = 1, cell%nlat
      th = cell%th(ilat)
      phi = cell%phi(ilon)

      cell%x(ilat,ilon,1) = r*sin(th)*cos(phi)
      cell%x(ilat,ilon,2) = r*sin(th)*sin(phi)
      cell%x(ilat,ilon,3) = r*cos(th)
    end do ! ilat
    end do ! ilon

    if (present(xc)) then
      do ii = 1, 3
        cell%x(:,:,ii) = cell%x(:,:,ii) + xc(ii)
      end do ! ii
    end if

  end subroutine RBC_MakeSphere
  
  ! Subroutine creatse an uncurved shape with the size of a sickle cell. This prolate spheroid shape is then
  ! run through a flow to create a curve in the shape. The resulting curved shape after the flow was then 
  ! exported and used as a sickle cell.
  ! Reference:
  !   Xiao Zhang, Wilbur A. Lam and Michael D. Graham
  !   Dynamics of deformable straight and curved prolate capsules in simple shear flow
  !   PHYSICAL REVIEW FLUIDS 5, 053101 (2020)
  subroutine RBC_MakeUncurvedSickleSpheroid(cell, r, xc)
  
    type(t_RBC) :: cell
    real(WP) :: r
    real(WP),optional :: xc(3)
    integer :: ilat, ilon, ii
    real(WP) :: th, phi
    real(WP) :: pol, eq, curv
  
    pol = 1.2 * r
    eq = 0.25 * r
    curv = 0.2 * r
  
    do ilon = 1, cell%nlon
    do ilat = 1, cell%nlat
      th = cell%th(ilat)
      phi = cell%phi(ilon)
      cell%x(ilat,ilon,1) = pol*sin(th)*cos(phi)
      cell%x(ilat,ilon,2) = eq*sin(th)*sin(phi)
      cell%x(ilat,ilon,3) = eq*cos(th)
  
      !curve capsule
      !cell%x(ilon, ilat, 2) = cell%x(ilon, ilat, 2) - (curv * (cell%x(ilon, ilat, 3)**2))
    end do ! ilat
    end do ! ilon
    if (present(xc)) then
      do ii = 1, 3
        cell%x(:,:,ii) = cell%x(:,:,ii) + xc(ii)
      end do ! ii
    end if
  
    end subroutine RBC_MakeUncurvedSickleSpheroid    

  ! Make Sickle Cell based on the shape given by Reference:
  !   Kviatkovsky, I. Zeidan, A. Yeheskely-Hayon, D. Shabad, E. L. Dann, E. J. & Yelin, D. 
  !   Measuring Sickle Cell morphology during blood flow
  !   Biomedical optics express, 8(3), 1996–2003
  !   doi: 10.1364/BOE.8.001996
  ! Unfortunately this shape doesn't work very well since points are ill-conditioned
  ! and diverge for this type of simulation, so we do not proceed with this sickle-cell model

  ! subroutine RBC_MakeSickle(cell, rad, xc)
  !   type(t_RBC) :: cell
  !   real(WP),optional :: xc(3)

  !   integer :: ilat, ilon, ii
  !   real(WP) :: th, phi, r_u, r_l, r, p, rad

  !   real(WP), dimension(0:5) :: a_l, a_u
  !   real(WP), dimension(2) :: b

  !   a_u = (/ 1.36 , -0.0403, 0.306, -0.00169, -0.0360, -0.0277 /)
  !   a_l = (/ -0.806, -0.1141, -0.00678, 0.00212, 0.0201,  0.0284 /)

  !   b = (/ 5.8, 3.05 /)
  !   p = 1.54

  !   do ilon = 1, cell%nlon
  !   do ilat = 1, cell%nlat
  !     th = cell%th(ilat)
  !     phi = cell%phi(ilon)

  !     r_u = RBC_SolveSickleRho(th, phi, a_u)
  !     r_l = RBC_SolveSickleRho(th, phi, a_l)
  !     r = min(r_u, r_l)
  !     r = r * rad / b(1)

  !     cell%x(ilat,ilon,1) = r*sin(th)*cos(phi)
  !     cell%x(ilat,ilon,2) = r*sin(th)*sin(phi)
  !     cell%x(ilat,ilon,3) = r*cos(th)

  !   end do ! ilat
  !   end do !ilon


  !   if (present(xc)) then
  !     do ii = 1, 3
  !       cell%x(:,:,ii) = cell%x(:,:,ii) + xc(ii)
  !     end do ! ii
  !   end if

  ! end subroutine RBC_MakeSickle

  ! Helper function for RBC_MakeSickle
  ! function RBC_SolveSickleRho(th, phi, a) result(r)
  !   real(WP) :: th, phi, r, ct2, st2, cp2, sp2
  !   real(WP),dimension(0:5) :: a ! parameter, coefficients for cartesian
  !   real(WP),dimension(5) :: quarticCoeffs
  !   complex(WP), dimension(4) :: all_roots
  !   integer :: i

  !   if (sin(th) .eq. 0) then
  !     r = a(0) / cos(th)
  !     return
  !   end if

  !   ct2 = cos(th) ** 2
  !   st2 = sin(th) ** 2
  !   cp2 = cos(phi) ** 2
  !   sp2 = sin(phi) ** 2

  !   quarticCoeffs(1) = a(0)
  !   quarticCoeffs(2) = -1 * cos(th)
  !   quarticCoeffs(3) = a(1) * cp2 * st2 + a(2) * sp2 * st2
  !   quarticCoeffs(4) = 0
  !   quarticCoeffs(5) = a(3) * ((cp2 * st2)**2) + a(4) * ((sp2 * st2)**2) + a(5) * (cp2 * sp2 *(st2 ** 2))

  !   call QuarticRoots(quarticCoeffs, all_roots)

  !   r = huge(r)
  !   do i = 1, 4
  !     if (REAL(AIMAG(all_roots(i))) .eq. 0 .and. REAL(REAL(all_roots(i))) .ge. 0) then
  !       r = min(r, REAL(REAL(all_roots(i))))
  !     end if
  !   end do

  ! end function RBC_SolveSickleRho
    
  
!**********************************************************************
! Make a spherical shaped cell for modeling a larger leukocyte --- r factor is hard coded here (for now at least)
! Arguments:
!  cell -- red blood cell
!  r -- spherical radius
!  xc -- center of the sphere

  subroutine RBC_MakePlatelet(cell, r, xc)
    type(t_RBC) :: cell
    real(WP) :: r
    real(WP),optional :: xc(3)
  
    real(WP) :: rplat

    rplat = 4.*3./(2.*2.82)/2.*r   ! Platelet expansion factor
    call RBC_MakeEllipsoid(cell, rplat, xc, 0.7) !0.39685)

  end subroutine RBC_MakePlatelet



!**********************************************************************
! Make a ellptically shaped cell
! Arguments:
!  cell -- red blood cell
!  r -- spherical radius
!  xc -- center of the sphere
  subroutine RBC_MakeEllipsoid(cell, r, xc, esn)
    type(t_RBC) :: cell
    real(WP) :: r,esn
    real(WP),optional :: xc(3)

    integer :: ilat, ilon, ii
    real(WP) :: th, phi

    do ilon = 1, cell%nlon
    do ilat = 1, cell%nlat
      th = cell%th(ilat)
      phi = cell%phi(ilon)

      cell%x(ilat,ilon,1) = 1./SQRT(esn)*r*sin(th)*cos(phi)
      cell%x(ilat,ilon,2) = 1./SQRT(esn)*r*sin(th)*sin(phi)
      cell%x(ilat,ilon,3) = esn*r*cos(th)
    end do ! ilat
    end do ! ilon

    if (present(xc)) then
      do ii = 1, 3
        cell%x(:,:,ii) = cell%x(:,:,ii) + xc(ii)
      end do ! ii
    end if

  end subroutine RBC_MakeEllipsoid
    
!**********************************************************************
! Make a bioconcave shaped red blood cell
! Arguments:
!  cell --
!  rc -- equivalent cell radius
!  xc -- center of the cell
! Reference:
!   -- C. Pozrikidis, Axisymmetric motion of a file of red blood cells through
!   capillaries, Physics of Fluids, 17, 031503 (2005).
  subroutine RBC_MakeBiConcave(cell, rc, xc)
    type(t_RBC) :: cell
    real(WP) :: rc
    real(WP),optional :: xc(3)

    integer :: ilat, ilon, ii
    real(WP) :: th, phi, r, z
    double precision :: alph = 1.3858189

    do ilat = 1, cell%nlat
      th = cell%th(ilat)

      z = rc * 0.5 * alph *(0.207 + 2.003*sin(th)**2 - 1.123*sin(th)**4) * cos(th)
      r = rc * alph * sin(th)

      do ilon = 1, cell%nlon
        phi = cell%phi(ilon)

        cell%x(ilat,ilon,1) = r*cos(phi)
        cell%x(ilat,ilon,2) = r*sin(phi)
        cell%x(ilat,ilon,3) = z
      end do ! ilon
    end do ! ilat

    if (present(xc)) then
      do ii = 1, 3
        cell%x(:,:,ii) = cell%x(:,:,ii) + xc(ii)
      end do ! ii
    end if

  end subroutine RBC_MakeBiConcave

!**********************************************************************
! Compute the geometry of a RBC
! Arguments:
!  cell -- 
  subroutine RBC_ComputeGeometry(cell)
    type(t_RBC) :: cell

    integer :: nlat, nlon, ilat, ilon
    real(WP),dimension(3) :: a1, a2, a3, a1_rcp, a2_rcp
    real(WP) :: a(2,2), a_rcp(2,2), detA, idetA, ds, b(2,2)
    real(WP),dimension(:,:,:),allocatable :: va, vb, a31, a32

    nlat = cell%nlat
    nlon = cell%nlon

    ! Allocate working arrays
    allocate(va(nlon/2+1,nlat,3), vb(nlon/2+1,nlat,3) )
    allocate(a31(nlat,nlon,3), a32(nlat,nlon,3) )

    call ShAnalGau(nlat, nlon, 3, cell%x, size(cell%x,1), size(cell%x,2), &
            va, vb, size(va,1), size(va,2), cell%wshags )
    call ShGradGau(nlat, nlon, 3, cell%a1, cell%a2, size(cell%a1,1), size(cell%a1,2),&
        va, vb, size(va,1), size(va,2), cell%wvhsgs )

    ! Surface tangential and normals
    do ilon = 1, nlon
    do ilat = 1, nlat
      a1 = cell%a1(ilat,ilon,:)
      a2 = cell%a2(ilat,ilon,:)

      a(1,1) = dot_product(a1, a1)
      a(1,2) = dot_product(a1, a2)
      a(2,1) = a(1,2)
      a(2,2) = dot_product(a2, a2)

      detA = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      iDetA = 1./detA

      a_rcp(1,1) =  iDeta*a(2,2)
      a_rcp(1,2) = -iDeta*a(1,2)
      a_rcp(2,1) = -iDeta*a(2,1)
      a_rcp(2,2) =  iDeta*a(1,1)

      a1_rcp = a_rcp(1,1)*a1 + a_rcp(1,2)*a2
      a2_rcp = a_rcp(2,1)*a1 + a_rcp(2,2)*a2

      ! Store results
      cell%a1_rcp(ilat,ilon,:) = a1_rcp
      cell%a2_rcp(ilat,ilon,:) = a2_rcp

      cell%a(ilat,ilon,:,:) = a
      cell%a_rcp(ilat,ilon,:,:) = a_rcp

      cell%detj(ilat,ilon) = sqrt(detA)/sin(cell%th(ilat))

      a3 = CrossProd(a1, a2)
      cell%a3(ilat,ilon,:) = a3/VecNorm(a3)
    end do ! ilat
    end do ! ilon

    ! Compute curvature tensor
    call ShAnalGau(nlat, nlon, 3, cell%a3, size(cell%a3,1), size(cell%a3,2), &
            va, vb, size(va,1), size(va,2), cell%wshags )
    call ShGradGau(nlat, nlon, 3, a31, a32, size(a31,1), size(a31,2), &
            va, vb, size(va,1), size(va,2), cell%wvhsgs )

    do ilon = 1, nlon
    do ilat = 1, nlat
      a1 = cell%a1(ilat,ilon,:)
      a2 = cell%a2(ilat,ilon,:)

      b(1,1) = -dot_product(a1, a31(ilat,ilon,:) )
      b(1,2) = -dot_product(a1, a32(ilat,ilon,:) )
      b(2,1) = -dot_product(a2, a31(ilat,ilon,:) )
      b(2,2) = -dot_product(a2, a32(ilat,ilon,:) )

      cell%b(ilat,ilon,:,:) = b 
    end do ! ilat
    end do ! ilon

    ! Compute centroid, area
    cell%xc = 0.
    cell%area = 0.

    do ilat = 1, nlat
    do ilon = 1, nlon
      ds = cell%detj(ilat,ilon)*cell%w(ilat)
      cell%area = cell%area + ds
      cell%xc = cell%xc + ds*cell%x(ilat,ilon,:)
    end do ! ilon
    end do ! ilat
    cell%xc = cell%xc/cell%area

    ! Compute volume
    cell%vol = 0.
    do ilon = 1, nlon
    do ilat = 1, nlat
      ds = cell%detj(ilat,ilon)*cell%w(ilat)
      cell%vol = cell%vol + ds*dot_product(cell%a3(ilat,ilon,:), cell%x(ilat,ilon,:)-cell%xc)
    end do ! ilat
    end do ! ilon
    cell%vol = THRD*cell%vol

    ! Mesh size
    cell%meshSize = sqrt(cell%area/cell%nlat**2)

    ! Deallocate working arrays
    deallocate(va, vb, a31, a32)

  end subroutine RBC_ComputeGeometry

!**********************************************************************
! Compute the covariant gradient of a vecotr field
! Arguments:
!  cell --
!  v(ilat,ilon,:) -- covariant vector component
!  dv(ilat,ilon,:,:) -- covariant gradient tensor
! Note:
!  -- v = v^i a_i
!  -- dv(i,j) = dv^i_j
!
!  -- The general formula to compute covariant derivatives does not apply 
!    here because of the discontinuity of the vector components at the two poles.
!
!  -- Instead, we transform the vector to a three-dimensional vector
!   in Eucledian space, compute the derivative of every Carteisan
!   component, and then project the derivative back onto the surface
  subroutine RBC_CovarGrad_Vec(cell, v, dv)
    type(t_RBC) :: cell
    real(WP) :: v(:,:,:), dv(:,:,:,:)

    integer :: nlat, nlon, ilat, ilon
    real(WP),dimension(3) :: a1, a2, a1_rcp, a2_rcp
    ! vc is the Cartesian representation of v in the 3D Eucledian space
    real(WP),dimension(:,:,:),allocatable :: vc, vc1, vc2, vca, vcb

    nlat = cell%nlat
    nlon = cell%nlon

    ! Allocate working arrays
    allocate(vc(nlat,nlon,3), vc1(nlat,nlon,3), vc2(nlat,nlon,3) )
    allocate(vca(nlon/2+1,nlat,3), vcb(nlon/2+1,nlat,3) )

    ! Transform the vector to representation in Eucledian space
    do ilat = 1, nlat
    do ilon = 1, nlon
      a1 = cell%a1(ilat,ilon,:)
      a2 = cell%a2(ilat,ilon,:)

      vc(ilat,ilon,:) = v(ilat,ilon,1)*a1 + v(ilat,ilon,2)*a2
    end do ! ilon
    end do ! ilat

    ! Compute the derivative of each Cartesian component
    call ShAnalGau(nlat, nlon, 3, vc, size(vc,1), size(vc,2), &
            vca, vcb, size(vca,1), size(vca,2), cell%wshags )
    call ShGradGau(nlat, nlon, 3, vc1, vc2, size(vc1,1), size(vc1,2),&
            vca, vcb, size(vca,1), size(vca,2), cell%wvhsgs )

    ! Project the derivative to the surface
    do ilat = 1, nlat
    do ilon = 1, nlon
      a1_rcp = cell%a1_rcp(ilat,ilon,:)
      a2_rcp = cell%a2_rcp(ilat,ilon,:)

      dv(ilat,ilon,1,1) = dot_product(a1_rcp, vc1(ilat,ilon,:) )
      dv(ilat,ilon,1,2) = dot_product(a1_rcp, vc2(ilat,ilon,:) )
      dv(ilat,ilon,2,1) = dot_product(a2_rcp, vc1(ilat,ilon,:) )
      dv(ilat,ilon,2,2) = dot_product(a2_rcp, vc2(ilat,ilon,:) )
    end do ! ilon
    end do ! ilat

    ! Deallocate working arrays
    deallocate(vc, vca, vcb, vc1, vc2)

  end subroutine RBC_CovarGrad_Vec

!**********************************************************************
! Compute the covariant gradient of a second-order tensor field
! Arguments:
!  cell --
!  t(ilat,ilon,:,:) -- a tensor field
!  dt(ilat,ilon,:,:,:) -- covariant gradient tensor
! Note:
!   -- t = t^{ij} e_i x e_j
!   -- t^{ij}_k
  subroutine RBC_CovarGrad_Tensor(cell, t, dt)
    type(t_RBC) :: cell
    real(WP) :: t(:,:,:,:), dt(:,:,:,:,:)

    integer :: nlat, nlon, ilat, ilon, i, j
    real(WP),dimension(:,:,:,:),allocatable :: tc, tc1, tc2, tca, tcb
    real(WP) :: a1(3), a2(3), a1_rcp(3), a2_rcp(3)
    real(WP) :: t_loc(2,2), t1_loc(2,2), t2_loc(2,2)
    real(WP) :: tc_loc(3,3), tc1_loc(3,3), tc2_loc(3,3)

    nlat = cell%nlat
    nlon = cell%nlon

    ! Allocate working arrays
    allocate(tc(nlat,nlon,3,3), tc1(nlat,nlon,3,3), tc2(nlat,nlon,3,3) )
    allocate(tca(nlon/2+1,nlat,3,3), tcb(nlon/2+1,nlat,3,3) )

    ! Transform to Eucledian space
    do ilat = 1, nlat
    do ilon = 1, nlon
      a1 = cell%a1(ilat,ilon,:)
      a2 = cell%a2(ilat,ilon,:)
      t_loc = t(ilat,ilon,:,:)
      tc_loc = 0.
      do i = 1, 3
      do j = 1, 3
    tc_loc(i,j) = t_loc(1,1)*a1(i)*a1(j) + t_loc(1,2)*a1(i)*a2(j) +&
              t_loc(2,1)*a2(i)*a1(j) + t_loc(2,2)*a2(i)*a2(j)
      end do ! i
      end do ! j

      tc(ilat,ilon,:,:) = tc_loc
    end do ! ilon
    end do ! ilat

    ! Compute derivative of each component
    call ShAnalGau(nlat, nlon, 9, tc, size(tc,1), size(tc,2), &
        tca, tcb, size(tca,1), size(tca,2), cell%wshags )
    call ShGradGau(nlat, nlon, 9, tc1, tc2, size(tc1,1), size(tc1,2), &
        tca, tcb, size(tca,1), size(tca,2), cell%wvhsgs )

    ! Project back to surface tensor form
    do ilat = 1, nlat
    do ilon = 1, nlon
      a1_rcp = cell%a1_rcp(ilat,ilon,:)
      a2_rcp = cell%a2_rcp(ilat,ilon,:)
      tc1_loc = tc1(ilat,ilon,:,:)
      tc2_loc = tc2(ilat,ilon,:,:)

      t1_loc(1,1) = dot_product(a1_rcp, matmul(tc1_loc, a1_rcp))
      t1_loc(1,2) = dot_product(a1_rcp, matmul(tc1_loc, a2_rcp))
      t1_loc(2,1) = dot_product(a2_rcp, matmul(tc1_loc, a1_rcp))
      t1_loc(2,2) = dot_product(a2_rcp, matmul(tc1_loc, a2_rcp))
      
      t2_loc(1,1) = dot_product(a1_rcp, matmul(tc2_loc, a1_rcp))
      t2_loc(1,2) = dot_product(a1_rcp, matmul(tc2_loc, a2_rcp))
      t2_loc(2,1) = dot_product(a2_rcp, matmul(tc2_loc, a1_rcp))
      t2_loc(2,2) = dot_product(a2_rcp, matmul(tc2_loc, a2_rcp))

      dt(ilat,ilon,:,:,1) = t1_loc
      dt(ilat,ilon,:,:,2) = t2_loc
    end do ! ilon
    end do ! ilat
    
    ! Deallocate working arrays
    deallocate(tc, tc1, tc2)
    deallocate(tca, tcb)

  end subroutine RBC_CovarGrad_Tensor

!**********************************************************************
! Project a quantity by a forward spherical harmonic transform
! followed by a backward transform
  subroutine RBC_SphProject(cell, nvar, u)
    type(t_rbc) :: cell
    integer :: nvar
    real(WP) :: u(cell%nlat,cell%nlon,nvar)

    integer :: nlat, nlon
    real(WP),allocatable :: va(:,:,:), vb(:,:,:)

    ! Allocate working arrays
    nlat = cell%nlat
    nlon = cell%nlon
    allocate(va(nlat,nlat,nvar), vb(nlat,nlat,nvar))

    ! Forward transfomr
    call ShAnalGau(nlat, nlon, nvar, u, size(u,1), size(u,2), &
            va, vb, size(va,1), size(va,2), cell%wshags)

    ! Filtering
    call ShFilter(nlat, nlon, nvar, va, vb, size(va,1), size(va,2), cell%nlat0, cell%nlon0)

    ! Backward transform
    call ShSynthGau(nlat, nlon, nvar, u, size(u,1), size(u,2), &
            va, vb, size(va,1), size(va,2), cell%wshsgs)

    ! Deallocate working arrays
    deallocate(va, vb)

  end subroutine RBC_SphProject

!**********************************************************************
! Compute the integral of a function over the cell
! Argument:
!  cell --
!  f -- values of the function at mesh points
  function RBC_Integral(cell, f) result(v)
    type(t_RBC) :: cell
    real(WP) :: f(:,:), v

    integer :: nlat, nlon, ilat, ilon
    real(WP) :: dv

    nlat = cell%nlat
    nlon = cell%nlon

    v = 0.
    do ilat = 1, nlat
      dv = 0.
      do ilon = 1, nlon
    dv = dv + f(ilat,ilon)*cell%detj(ilat,ilon)
      end do ! ilon

      v = v + dv*cell%w(ilat)
    end do ! ilat

  end function RBC_Integral

!**********************************************************************
  subroutine Rbc_BuildSurfaceSource(cell, xFlag, fFlag, gFlag)
    type(t_Rbc) :: cell
    logical,optional :: xFlag, fFlag, gFlag

    integer :: nlat0, nlon0, nlat, nlon
    real(WP),dimension(:,:,:),allocatable :: v, va, vb
    real(WP),dimension(:,:,:),allocatable :: FdetJ, GdetJ
    integer :: ii

    nlat0 = cell%nlat0
    nlon0 = cell%nlon0

    nlat = cell%nlat
    nlon = cell%nlon

    ! Allocate working arrays
    allocate(v(0:nlat,nlon,3) )
    allocate(va(nlat+1,nlat+1,3), vb(nlat+1,nlat+1,3) )

    ! Build spline representation for x, a3, and detj
    if (present(xFlag)) then
      if (xFlag) then
    ! x
        call ShAnalGau(nlat, nlon, 3, cell%x, size(cell%x,1), size(cell%x,2), &
        va, vb, size(va,1), size(va,2), cell%wshags )
    call ShFilter(nlat, nlon, 3, va, vb, size(va,1), size(va,2), nlat0, nlon0)
        call ShSynthEqu(nlat+1, nlon, 3, v, size(v,1), size(v,2), &
        va, vb, size(va,1), size(va,2), cell%wshses )
        call Spline_Build_On_Sphere(cell%spln_x, v)


    ! a3
    call ShAnalGau(nlat, nlon , 3, cell%a3, size(cell%a3,1), size(cell%a3,2), &
        va, vb, size(va,1), size(va,2), cell%wshags )
        call ShFilter(nlat, nlon, 3, va, vb, size(va,1), size(va,2), nlat0, nlon0 )
        call ShSynthEqu(nlat+1, nlon, 3, v, size(v,1), size(v,2), &
        va, vb, size(va,1), size(va,2), cell%wshses )
        call Spline_Build_On_Sphere(cell%spln_a3, v)

    ! detJ
    call ShAnalGau(nlat, nlon, 1, cell%detj, size(cell%detj,1), size(cell%detj,2), &
        va, vb, size(va,1), size(va,2), cell%wshags )
        call ShFilter(nlat, nlon, 1, va, vb, size(va,1), size(va,2), nlat0, nlon0 )
        call ShSynthEqu(nlat+1, nlon, 1, v, size(v,1), size(v,2), &
        va, vb, size(va,1), size(va,2), cell%wshses )
        call Spline_Build_On_Sphere(cell%spln_detJ, v(:,:,1))
      end if
    end if

    ! f*detJ
    if (present(fFlag)) then
      if (fFlag) then
    ! Spline representation
    allocate(FdetJ(nlat,nlon,3) )
    do ii = 1, 3
      FdetJ(:,:,ii) = cell%f(:,:,ii) * cell%detj
        end do ! ii

        call ShAnalGau(nlat, nlon, 3, FdetJ, size(FdetJ,1), size(FdetJ,2), &
        va, vb, size(va,1), size(va,2), cell%wshags )
        call ShFilter(nlat, nlon, 3, va, vb, size(va,1), size(va,2), nlat0, nlon0 )
        call ShSynthEqu(nlat+1, nlon, 3, v, size(v,1), size(v,2), &
        va, vb, size(va,1), size(va,2), cell%wshses )
        call Spline_Build_On_Sphere(cell%spln_FdetJ, v)

    deallocate(FdetJ)
      end if
    end if

    ! g*detJ
    if (present(gFlag)) then
      if (gFlag) then
    ! Spline representation
    allocate(GdetJ(nlat,nlon,3) )
    do ii = 1, 3
      gdetj(:,:,ii) = cell%g(:,:,ii) * cell%detj
        end do ! ii

        call ShAnalGau(nlat, nlon, 3, GdetJ, size(GdetJ,1), size(GdetJ,2), &
        va, vb, size(va,1), size(va,2), cell%wshags )
        call ShFilter(nlat, nlon, 3, va, vb, size(va,1), size(va,2), nlat0, nlon0 )
        call ShSynthEqu(nlat+1, nlon, 3, v, size(v,1), size(v,2), &
        va, vb, size(va,1), size(va,2), cell%wshses )
        call Spline_Build_On_Sphere(cell%spln_GdetJ, v)

    deallocate(GdetJ)
      end if
    end if

    ! Deallocate working arrays
    deallocate(v, va, vb)

  end subroutine Rbc_BuildSurfaceSource

!**********************************************************************
! Compute the force needed for force balance on the cell
! Arguments:
!  cell -- deformed cell
!  cellRef -- reference cell
!  f -- residual force, in Cartesian coordinates
  subroutine Shell_ResForce(cell, cellRef, f)
    type(t_rbc) :: cell, cellRef
    real(WP) :: f(:,:,:)

    integer :: nlat, nlon, ilat, ilon
    real(WP),allocatable :: bnd(:,:,:,:), dbnd(:,:,:,:,:), &
        tau_t(:,:,:,:), dtau_t(:,:,:,:,:), tau_n(:,:,:), dtau_n(:,:,:,:)
    real(WP) :: a_rcp(2,2), b(2,2), f_tmp(3)
    real(WP) :: fnTot

    ! Initialize
    f = 0.

    ! Allocate temoprary arrays
    nlat = cell%nlat
    nlon = cell%nlon
    allocate(bnd(nlat,nlon,2,2), dbnd(nlat,nlon,2,2,2) )
    allocate(tau_t(nlat,nlon,2,2), dtau_t(nlat,nlon,2,2,2) )
    allocate(tau_n(nlat,nlon,2), dtau_n(nlat,nlon,2,2) )

    ! Compute in-plane elastic stress and bending moment
    call Shell_ElasStrs(cell, cellRef, tau_t)
    call Shell_Bend(cell, cellRef, bnd)

    ! Compute the asymmetric part of the in-plane stress
    ! To be added

    ! Compute normal stress
    call RBC_CovarGrad_Tensor(cell, bnd, dbnd)
    do ilat = 1, nlat
    do ilon = 1, nlon
      tau_n(ilat,ilon,1) = dbnd(ilat,ilon,1,1,1) + dbnd(ilat,ilon,2,1,2)
      tau_n(ilat,ilon,2) = dbnd(ilat,ilon,1,2,1) + dbnd(ilat,ilon,2,2,2)
    end do ! ilon
    end do ! ilat

    ! Compute tangential residual force
    call RBC_CovarGrad_Tensor(cell, tau_t, dtau_t)
    do ilat = 1, nlat
    do ilon = 1, nlon
      f_tmp(1) = - dtau_t(ilat,ilon,1,1,1) - dtau_t(ilat,ilon,2,1,2)
      f_tmp(2) = - dtau_t(ilat,ilon,1,2,1) - dtau_t(ilat,ilon,2,2,2)

      a_rcp = cell%a_rcp(ilat,ilon,:,:)
      b = cell%b(ilat,ilon,:,:)
      f_tmp(1:2) = f_tmp(1:2) + matmul(a_rcp, matmul(b, tau_n(ilat,ilon,:)))

      f(ilat,ilon,1:2) = f_tmp(1:2)
    end do ! ilon
    end do ! ilat

    ! Compute the normal residual force
    call RBC_CovarGrad_Vec(cell, tau_n, dtau_n)
    do ilat = 1, nlat
    do ilon = 1, nlon
      b = cell%b(ilat,ilon,:,:)

      f(ilat,ilon,3) = -dtau_n(ilat,ilon,1,1) - dtau_n(ilat,ilon,2,2) &
      - sum(tau_t(ilat,ilon,:,:)*b)
    end do ! ilon
    end do ! ilat

    ! Subtract the average normal force
    fnTot = 0.
    do ilon = 1, nlon
    do ilat = 1, nlat
      fnTot = fnTot + f(ilat,ilon,3) * cell%detJ(ilat,ilon) * cell%w(ilat)
    end do ! ilat
    end do ! ilon
    f(:,:,3) = f(:,:,3) - fnTot/cell%area

    ! Convert the force to Cartesian coordinate
    do ilat = 1, nlat
    do ilon = 1, nlon
      f(ilat,ilon,:) = f(ilat,ilon,1)*cell%a1(ilat,ilon,:) + &
      f(ilat,ilon,2)*cell%a2(ilat,ilon,:) + &
      f(ilat,ilon,3)*cell%a3(ilat,ilon,:)
    end do ! ilon
    end do ! ilat

    ! Deallocate working arrays
    deallocate(bnd, dbnd )
    deallocate(tau_t, dtau_t )
    deallocate(tau_n, dtau_n )

  end subroutine Shell_ResForce

!**********************************************************************
! Compute the elastic stress due to stretches
! Arguments:
!  cell -- cell in current configurations
!  cellRef -- cell in a reference state
!  strs -- elastic stress tensor in contravariant form
  subroutine Shell_ElasStrs(cell, cellRef, tau)
    type(t_rbc) :: cellRef, cell
    real(WP) :: tau(:,:,:,:)

    integer :: nlat, nlon, ilat, ilon, i, j
    real(WP) :: V2(2,2), JS, lbd1, lbd2, tau_tmp(2,2)
    real(WP) :: ES, ED

    ES = cell%ES
    ED = cell%ED

    nlat = cell%nlat
    nlon = cell%nlon

    do ilat = 1, nlat
    do ilon = 1, nlon
      ! V2 is the strech tensor
      V2 = matmul(cell%a(ilat,ilon,:,:), cellRef%a_rcp(ilat,ilon,:,:) )

      lbd1 = V2(1,1) + V2(2,2) - 2.0
      lbd2 = V2(1,1)*V2(2,2) - V2(1,2)*V2(2,1) - 1.0
      JS = sqrt(V2(1,1)*V2(2,2) - V2(1,2)*V2(2,1))

      ! Take the covariant form of V2
      V2 = cellRef%a_rcp(ilat,ilon,:,:)

      tau_tmp = 0.5*ES/JS*(lbd1 + 1.0)*V2 &
              + 0.5*JS*(ED*lbd2 - ES)*cell%a_rcp(ilat,ilon,:,:)

      tau(ilat,ilon,:,:) = tau_tmp
    end do ! ilon
    end do ! ilat

  end subroutine Shell_ElasStrs

!**********************************************************************
! Compute the in-surface bending momentum
! Arguments:
! cell - deformed cell
! cellRef -- reference cell
! bnd -- bending moment tensor, covariant coefficients
  subroutine Shell_Bend(cell, cellRef, bnd)
    type(t_rbc) :: cell, cellRef
    real(WP) :: bnd(:,:,:,:)

    integer :: nlat, nlon, ilat, ilon
    real(WP),dimension(2,2) :: b, bRef

    nlat = cell%nlat
    nlon = cell%nlon

    do ilat = 1, nlat
    do ilon = 1, nlon
      b = matmul(cell%a_rcp(ilat,ilon,:,:), cell%b(ilat,ilon,:,:))
      bref = matmul(cellRef%a_rcp(ilat,ilon,:,:), cellRef%b(ilat,ilon,:,:))

      bnd(ilat,ilon,:,:) = -cell%EB*matmul(b - bref, cell%a_rcp(ilat,ilon,:,:))
    end do ! ilon
    end do ! ilat

  end subroutine Shell_Bend

!**********************************************************************

end module ModRbc
