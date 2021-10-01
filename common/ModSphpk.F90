! Wrapper for spherepack library
module ModSphpk

  use ModDataTypes

  implicit none

  public ::  ShAnalGau_Init, &
    ShAnalGau, &
    ShAnalEqu_Init, &
    ShAnalEqu, &
    ShSynthGau_Init, &
    ShSynthGau, &
    ShSynthEqu_Init, &
    ShSynthEqu, &
    ShGradGau_Init, &
    ShGradGau, &
    ShFilter

contains

!**********************************************************************
! Compute the working array for scalar spherical harmonic analysis
! Arguments:
!   nlat, nlon --
!   wshags -- the working array that can be reused
  subroutine ShAnalGau_Init(nlat, nlon, wshags)
    integer :: nlat, nlon
    real(WP),pointer :: wshags(:)

    integer :: l1, l2, lshags
    integer :: lwork, ldwork
    real(WP),allocatable ::  work(:), dwork(:)
    character(*),parameter :: func_name = 'ShAnalGau_Init'
    integer :: ierr

    ! Allocate working arrays
    l1 = min(nlat, nlon/2+1)
    l2 = (nlat + 1)/2
    lshags = nlat*(3*(l1+l2)-2) + (l1-1)*(l2*(2*nlat-l1)-3*l1)/2 + nlon + 15
    allocate(wshags(lshags) )

    lwork = 4*nlat*(nlat + 2) + 2
    allocate(work(lwork) )

    ldwork = nlat*(nlat+4)
    allocate(dwork(ldwork) )

    call shagsi(nlat, nlon, wshags, lshags, work, lwork, dwork, ldwork, ierr)

    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, nlat,nlon
      print *, 'Error calling shagsi, ierror = ', ierr
      stop
    end if

    ! Deallocate working arrays
    deallocate(work, dwork)
    
  end subroutine ShAnalGau_Init


!**********************************************************************
! Spherical harmonic decomposition on a Gaussian mesh
! Arguments:
!  nlat, nlon -- number of latitudinal, longitudinal points
!  nt -- number of variables
!  g(:,:,:) -- input array
!  idg, jdg -- 1st and 2nd dimensions of g (latitude and longitude)
!  a, b -- spherical harmonic coefficients
!  ida, jda -- 1st and 2nd dimensions of a
!  wshags -- working array that must be initialized by calling
!        ShAnalGau_Init
  subroutine ShAnalGau(nlat, nlon, nt, g, idg, jdg, &
        a, b, ida, jda, wshags)
    integer :: nlat, nlon, nt
    real(WP),dimension(idg,jdg,1) :: g
    integer :: idg, jdg
    real(WP),dimension(ida,jda,1) :: a, b
    integer :: ida, jda
    real(WP) :: wshags(:)

    integer,parameter :: isym = 0
    real(WP),allocatable ::  work(:)
    integer,save :: lwork
    integer :: ierr
    character(*),parameter :: func_name = 'ShAnalGau'

    ! Allocate working arrays
    lwork = nlat*nlon*(nt + 1)
    allocate(work(lwork))

    call shags(nlat, nlon, isym, nt, g, size(g,1), size(g,2), &
      a, b, size(a,1), size(a,2), wshags, size(wshags), work, lwork, ierr)

    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, 'Error calling shags, ierror = ', ierr
      stop
    end if

    ! Deallocate working arrays
    deallocate(work)

  end subroutine ShAnalGau


!**********************************************************************
! Compute the working array for scalar spherical harmonic analysis
! on a uniform mesh
! Arguments:
!   nlat, nlon --
!   wshags -- the working array that can be reused
  subroutine ShAnalEqu_Init(nlat, nlon, wshaes)
    integer :: nlat, nlon
    real(WP),pointer :: wshaes(:)

    integer :: l1, l2, lshaes
    integer :: lwork, ldwork
    real(WP),allocatable ::  work(:), dwork(:)
    character(*),parameter :: func_name = 'ShAnalEqu_Init'
    integer :: ierr

    ! Allocate working arrays
    l1 = min(nlat, (nlon+2)/2)
    l2 = (nlat + 1)/2
    lshaes = (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
    allocate(wshaes(lshaes) )

    l1 = min(nlat, (nlon+2)/2)
    l2 = (nlat + 1)/2
    lwork = 5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
    allocate(work(lwork) )

    ldwork = nlat + 1
    allocate(dwork(ldwork) )

    call shaesi(nlat, nlon, wshaes, lshaes, work, lwork, dwork, ldwork, ierr)

    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, 'Error calling shaesi, ierror = ', ierr
      stop
    end if

    ! Deallocate working arrays
    deallocate(work, dwork)
    
  end subroutine ShAnalEqu_Init

!**********************************************************************
! Spherical harmonic decomposition on a equally spaced mesh
! Arguments:
!  nlat, nlon -- number of latitudinal, longitudinal points
!  nt -- number of variables
!  g(:,:,:) -- input array
!  idg, jdg -- 1st and 2nd dimensions of g (latitude and longitude)
!  a, b -- spherical harmonic coefficients
!  ida, jda -- 1st and 2nd dimensions of a
!  wshags -- working array that must be initialized by calling
!        ShAnalGau_Init
!
  subroutine ShAnalEqu(nlat, nlon, nt, g, idg, jdg, &
        a, b, ida, jda, wshaes)
    integer :: nlat, nlon, nt
    real(WP),dimension(idg,jdg,1) :: g
    integer :: idg, jdg
    real(WP),dimension(ida,jda,1) :: a, b
    integer :: ida, jda
    real(WP) :: wshaes(:)

    integer,parameter :: isym = 0
    real(WP),allocatable ::  work(:)
    integer,save :: lwork
    integer :: ierr
    character(*),parameter :: func_name = 'ShAnalEqu'

    ! Allocate working arrays
    lwork = (nt+1)*nlat*nlon
    allocate(work(lwork))

    call shaes(nlat, nlon, isym, nt, g, size(g,1), size(g,2), &
      a, b, size(a,1), size(a,2), wshaes, size(wshaes), work, lwork, ierr)

    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, 'Error calling shaes, ierror = ', ierr
      stop
    end if

    ! Deallocate working arrays
    deallocate(work)

  end subroutine ShAnalEqu

!**********************************************************************
  subroutine ShSynthGau_Init(nlat, nlon, wshsgs)
    integer :: nlat, nlon
    real(WP),pointer :: wshsgs(:)

    integer :: l1, l2, lshsgs
    integer :: lwork, ldwork
    real(WP),allocatable :: work(:), dwork(:)
    character(*),parameter :: func_name = 'ShSynthGau_Init'
    integer :: ierr

    ! Allocate working arrays
    l1 = min(nlat, nlon/2+1)
    l2 = (nlat + 1)/2
    lshsgs = nlat*(3*(l1+l2)-2) + (l1-1)*(l2*(2*nlat-l1)-3*l1)/2 + nlon + 15
    allocate(wshsgs(lshsgs) )

    lwork = 4*nlat*(nlat+2)+2
    allocate(work(lwork) )

    ldwork = nlat*(nlat+4)
    allocate(dwork(ldwork) )

    call shsgsi(nlat, nlon, wshsgs, lshsgs, work, lwork, dwork, ldwork,ierr)
    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, 'Error when calling shsgsi, ierror = ', ierr
      stop
    end if

    ! Deallocate working arrays
    deallocate(work, dwork)

  end subroutine ShSynthGau_Init

!**********************************************************************
! Systhesis 
! Arguments:
!  nlat, nlon, nt --
!  g -- physical variables
!  idg, jdg -- 1st and 2nd dimensions of g
!  a, b -- spherical harmonic coefficients
!  ida, idb -- 1st and 2nd dimensions of a
! Note:
!  -- See document for spherepack subroutines shsgs and shsgi
  subroutine ShSynthGau(nlat, nlon, nt, g, idg, jdg, a, b, ida, jda, wshsgs)
    integer :: nlat, nlon, nt
    real(WP),dimension(idg,jdg,1) :: g
    integer :: idg, jdg
    real(WP),dimension(ida,jda,1) :: a, b
    integer :: ida, jda
    real(WP) :: wshsgs(:)

    integer,parameter :: isym = 0
    integer :: lwork
    real(WP),allocatable:: work(:)
    integer :: ierr
    character(*),parameter :: func_name = 'ShSynthGau'

    ! Allocate working arrays
    lwork = nlat*nlon*(nt + 1)
    allocate(work(lwork) )

    call shsgs(nlat, nlon, isym, nt, g, idg, jdg, &
      a, b, ida, jda, wshsgs, size(wshsgs), work, lwork, ierr)

    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, 'Error calling shsgs, ierror = ', ierr
      stop
    end if

    ! Deallocate working arrays
    deallocate(work)

  end subroutine ShSynthGau

!**********************************************************************
  subroutine ShSynthEqu_Init(nlat, nlon, wshses)
    integer :: nlat, nlon
    real(WP),pointer :: wshses(:)

    integer :: l1, l2, lshses, lwork, ldwork
    real(WP),allocatable :: work(:), dwork(:)
    integer :: ierr
    character(*),parameter :: func_name = 'ShSynthEqu_Init'
    
    ! Allocate working arrays
    l1 = min(nlat, nlon/2+1)
    l2 = (nlat + 1)/2
    lshses = (l1*l2*(nlat+nlat-l1+1))/2+nlon+15
    allocate(wshses(lshses))

    lwork = 5*nlat*l2+3*((l1-2)*(nlat+nlat-l1-1))/2
    allocate(work(lwork))

    ldwork = nlat + 1
    allocate(dwork(ldwork) )

    call shsesi(nlat, nlon, wshses, lshses, work, lwork, dwork, ldwork, ierr)

    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, 'Error calling shsesi, ierror = ', ierr
      stop
    end if

    ! Deallocate working arrays
    deallocate(work, dwork)

  end subroutine ShSynthEqu_Init

!**********************************************************************
  subroutine ShSynthEqu(nlat, nlon, nt, g, idg, jdg, a, b, ida, jda, wshses)
    integer :: nlat, nlon, nt
    real(WP),dimension(idg,jdg,1) :: g
    integer :: idg, jdg
    real(WP),dimension(ida,jda,1) :: a, b
    integer :: ida, jda
    real(WP) :: wshses(:)

    integer,parameter :: isym = 0
    integer :: lwork
    real(WP),allocatable:: work(:)
    character(*),parameter :: func_name = 'ShSynthEqu'
    integer :: ierr

    ! Allocate working arrays
    lwork = nlat*nlon*(nt + 1)
    allocate(work(lwork) )

    call shses(nlat, nlon, isym, nt, g, idg, jdg, &
      a, b, ida, jda, wshses, size(wshses), work, lwork, ierr)

    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, 'Error calling shses, ierror = ', ierr
      stop
    end if

    ! Deallocate working arrays
    deallocate(work)

  end subroutine ShSynthEqu

!**********************************************************************
  subroutine ShGradGau_Init(nlat, nlon, wvhsgs)
    integer :: nlat, nlon
    real(WP),pointer :: wvhsgs(:)

    integer :: imid, lmn, lvhsgs, ldwork
    real(WP),allocatable :: dwork(:)
    integer :: ierr
    character(*),parameter :: func_name = 'ShGradGau_Init'

    ! Allocate working arrays
    imid = (nlat+1)/2
    lmn = (nlat*(nlat+1))/2
    lvhsgs = 2*(imid*lmn+nlat) + nlon + 15 
    allocate(wvhsgs(lvhsgs))

    ldwork = nlat*(3*nlat+8)+1
    allocate(dwork(ldwork) )

    call vhsgsi(nlat, nlon, wvhsgs, lvhsgs, dwork, ldwork, ierr)
    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, 'Error calling vhsgsi, ierror = ', ierr
      stop
    end if

    ! Deallocate working arrays
    deallocate(dwork)

  end subroutine ShGradGau_Init

!**********************************************************************
! Compute the gradient on a sphere
! Arguments:
!  v -- d g/d th, th is colatitude angle
!  w -- d g/d phi, w is longitude angle
!  a, b -- coefficients of spherical harmonic expansions of a scalar field
  subroutine ShGradGau(nlat, nlon, nt, v, w, idv, jdv, a, b, ida, jda, wvhsgs)
    integer :: nlat, nlon, nt
    real(WP),dimension(idv,jdv,1) :: v, w
    integer :: idv, jdv
    real(WP),dimension(ida,jda,1) :: a, b
    integer :: ida, jda
    real(WP) :: wvhsgs(:)
    
    integer,parameter :: isym = 0
    integer :: l1, l2, lwork
    real(WP),allocatable :: work(:)
    real(WP) :: th(nlat), weight(nlat)
    integer :: ilat
    integer :: ierr
    character(*),parameter :: func_name = 'ShGradGau'

    ! Allocate working arrays
    l1 = min(nlat, (nlon+1)/2)
    l2 = (nlat+1)/2
    lwork = nlat*((2*nt+1)*nlon+2*l1*nt+1)
    allocate(work(lwork) )

    call gradgs(nlat, nlon, isym, nt, v, w, idv, jdv, &
      a, b, ida, jda, wvhsgs, size(wvhsgs), work, lwork, ierr)

    if (ierr /= 0) then
      print *, 'Subroutine ', func_name
      print *, 'Error calling gradgs, ierror = ', ierr
      stop
    end if

    ! Add back the sin(th) factor
    call gaqd(nlat, th, weight, work, lwork, ierr)
    do ilat = 1, nlat
      w(ilat,:,1:nt) = sin(th(ilat)) * w(ilat,:,1:nt)
    end do ! ilat

    ! Deallocate working arrays
    deallocate(work)

  end subroutine ShGradGau

!**********************************************************************
! Filtering
! Note:
!  -- The filtering is such that a(i,j) /= 0 only for 1 <= i <= nlat0
  subroutine ShFilter(nlat, nlon, nt, a, b, ida, jda, nlat0, nlon0)
    integer :: nlat, nlon, nt
    real(WP),dimension(ida,jda,1) :: a, b
    integer :: ida, jda
    integer :: nlat0, nlon0

    a(:,nlat0+1:,1:nt) = 0.
    b(:,nlat0+1:,1:nt) = 0.

  end subroutine ShFilter

!**********************************************************************

end module ModSphpk
