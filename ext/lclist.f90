!----------------------------------------------------------------------
!      lclist.f90 - implementation of a linked cell neighbor list
!----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2016 Nongnuch Artrith and Alexander Urban
!+
!+ This program is free software: you can redistribute it and/or modify
!+ it under the terms of the GNU General Public License as published by
!+ the Free Software Foundation, either version 3 of the License, or
!+ (at your option) any later version.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!+ General Public License for more details.
!+
!+ You should have received a copy of the GNU General Public License
!+ along with this program.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------

module lclist

  !--------------------------------------------------------------------!
  ! A simple, but universal neighbourlist implementation without any   !
  ! restrictions on the lattice vectors or on the number of atoms in   !
  ! the system.  It is equally well suited for small unit cells, where !
  ! |lattice vectors| << R_cut, and for unit cells that are much       !
  ! larger than the interaction radius.                                !
  !                                                                    !
  ! If the simulation cell is large wrt. the interaction cut-off, it   !
  ! is partitioned into cells and the information about the atoms in   !
  ! each cell is stored in a linked cell list.                         !
  !--------------------------------------------------------------------!
  ! 2011-11-05 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use sortlib, only: argsort

  implicit none
  save

  public  :: lcl_init,               &
             lcl_final,              &
             lcl_print_info,         &
             lcl_nmax_cell,          &
             lcl_nmax_nblist,        &
             lcl_nmax_nbdist,        &
             lcl_nblist,             &
             lcl_nbdist,             &
             lcl_nbdist_cart

  private :: cell_assign_atoms,      &
             cell_multiples,         &
             cell_get,               &
             cell_add,               &
             wrap_cooLatt,           &
             cell_volume,            &
             inner_sphere,           &
             vproduct,               &
             translation_vectors


  double precision, parameter, private :: PI = 3.1415926535897931d0

  !--------------------------------------------------------------------!

  logical,                                         private :: pbc
  double precision, dimension(3,3),                private :: latticeVec
  double precision, dimension(3,3),                private :: gridVec
  integer,                                         private :: nAtoms
  integer,          dimension(:),   pointer,       private :: atomType
  double precision, dimension(:,:), pointer,       private :: cooLatt

  double precision,                                private :: r_min
  double precision,                                private :: r_max
  integer,          dimension(3),                  private :: nCells
  integer,          dimension(:,:,:), allocatable, private :: cell
  integer,          dimension(:),     allocatable, private :: atomList
  integer,          dimension(:,:),   allocatable, private :: cellList

  integer,                                         private :: nCvecs
  integer,          dimension(:,:),   allocatable, private :: Cvec

  integer,                                         private :: nTvecs
  integer,          dimension(:,:),   allocatable, private :: Tvec
  double precision, dimension(:),     allocatable, private :: TvecLen

  integer,                                         private :: nmax_cell
  integer,                                         private :: nmax_nblist
  integer,                                         private :: nmax_nbdist

  logical,                                         private :: isInit = .false.

contains

  subroutine lcl_init(r_min_in, r_max_in, latticeVec_in, &
                      nAtoms_in, atomType_in, cooLatt_in, pbc_in)

    implicit none

    double precision,                                 intent(in) :: r_min_in
    double precision,                                 intent(in) :: r_max_in
    double precision, dimension(3,3),                 intent(in) :: latticeVec_in
    integer,                                          intent(in) :: nAtoms_in
    integer,          dimension(nAtoms_in),   target, intent(in) :: atomType_in
    double precision, dimension(3,nAtoms_in), target, intent(in) :: cooLatt_in
    logical, optional,                                intent(in) :: pbc_in

    if (isInit) then
       write(0,*) "Error: module already initialized in `lclist_init'."
       return
    end if

    r_min           = r_min_in
    r_max           = r_max_in
    nAtoms          = nAtoms_in
    latticeVec(:,:) = latticeVec_in(:,:)
    atomType        => atomType_in(:)
    cooLatt         => cooLatt_in(:,:)
    if (present(pbc_in)) then
       pbc = pbc_in
    else
       pbc = .true.
    end if

    ! get optimal divison of the unit cell:
    call cell_multiples(r_max, latticeVec, gridVec, nCells)
    allocate(cell(nCells(3),nCells(2),nCells(1)), atomList(nAtoms), &
             cellList(3,nAtoms))

    ! assign atoms to cells:
    call cell_assign_atoms(nCells, nAtoms, cooLatt, cell, atomList, cellList)

    ! set up half start of vectors pointing to cells within range of r_max:
    nCvecs = 0
    call translation_vectors(r_max, gridVec, nCvecs, Cvec, nc=nCells)
    allocate(Cvec(3,nCvecs))
    call translation_vectors(r_max, gridVec, nCVecs, cVec, nc=nCells)

    ! set up half start of translation vectors pointing to periodic
    ! images of the unit cell within range of r_max:
    nTvecs = 0
    if (pbc) then
       call translation_vectors(r_max, latticeVec, nTvecs, Tvec, Tnorm=TvecLen)
       allocate(Tvec(3,nTvecs), TvecLen(nTvecs))
       call translation_vectors(r_max, latticeVec, nTvecs, Tvec, Tnorm=TvecLen)
    end if

    isInit = .true.

    ! max. number of neighbours and max. number of atoms per cell:
    nmax_cell   = lcl_nmax_cell()
    nmax_nblist = lcl_nmax_nblist()
    nmax_nbdist = lcl_nmax_nbdist(r_min, r_max)

  end subroutine lcl_init

  !--------------------------------------------------------------------!

  subroutine lcl_final()

    implicit none

    if (.not. isInit) return

    if (allocated(cell)) deallocate(cell, atomList, cellList)
    if (allocated(Cvec)) deallocate(Cvec)
    if (allocated(Tvec)) deallocate(Tvec, TvecLen)

    cooLatt  => null()
    atomType => null()

    isInit = .false.

  end subroutine lcl_final

  !--------------------------------------------------------------------!
  !       print information about the linked cell list to stdout       !
  !--------------------------------------------------------------------!

  subroutine lcl_print_info()

    implicit none

    write(*,*) 'Linked Cell List (Neighbourlist)'
    write(*,*) '--------------------------------'
    write(*,*)

    if (.not. isInit) then
       write(*,*) "Module not initialized."
       write(*,*)
       return
    end if

    write(*,'(1x,"Number of atoms             : ",I10)')      nAtoms
    write(*,'(1x,"Minimal distance            : ",ES10.3)')   r_min
    write(*,'(1x,"Cut-off radius              : ",ES10.3)')   r_max
    write(*,'(1x,"Number of cells             : ",3(I4,2x))') nCells
    write(*,'(1x,"Cells within cut-off        : ",I10)')      2*nCvecs + 1
    write(*,'(1x,"Translation vectors         : ",I10)')      2*nTvecs + 1
    write(*,'(1x,"Max. atoms per cell         : ",I10)')      nmax_cell
    write(*,'(1x,"Max. possible neighbour IDs : ",I10)')      nmax_nblist
    write(*,'(1x,"Max. real neighbours        : ",I10)')      nmax_nbdist
    write(*,*)

  end subroutine lcl_print_info

  !--------------------------------------------------------------------!
  !        max. number of neighbour candidates within cut-off          !
  !--------------------------------------------------------------------!

  function lcl_nmax_nblist() result(nmax)

    implicit none

    integer          :: nmax

    nmax = min(lcl_nmax_cell()*(2*nCvecs + 1), nAtoms-1)

  end function lcl_nmax_nblist

  !--------------------------------------------------------------------!
  !          max. number of (real) neighbours within cut-off           !
  !--------------------------------------------------------------------!

  function lcl_nmax_nbdist(rmin, rmax) result(nmax)

    implicit none

    double precision, intent(in) :: rmin, rmax
    integer                      :: nmax

    double precision :: V_atom, V_cut

    V_atom = (0.5d0*rmin)**3
    V_cut  = (rmax+0.5d0*rmin)**3

    ! max number of atoms assuming close packing
    ! pi/(s*sqrt(2)) ~ 0.7405
    nmax = ceiling(V_cut/V_atom*0.7405d0)

  end function lcl_nmax_nbdist

  !--------------------------------------------------------------------!
  !           estimated max. average number of atoms per cell          !
  !--------------------------------------------------------------------!

  function lcl_nmax_cell() result(nmax)

    implicit none

    integer          :: nmax
    double precision :: V_cell, V_atom

    if (.not. isInit) then
       write(0,*) "Error: module not initialized in `max_atoms_per_cell'."
       stop
    end if

    V_cell = cell_volume(gridVec)
    ! effective atom volume in a densly packed crystal structure:
    ! 3*sqrt(2)/pi * 4/3*pi * r^3 = 4*sqrt(2) * r^3
    V_atom = sqrt(32.0d0)*(0.5d0*r_min)**3
    ! wrong?! assuming the also partial spheres in the cell.
    ! V_atom = 8.0d0/27.0d0*(0.5d0*r_min)**3
    nmax   = ceiling(V_cell/V_atom)

  end function lcl_nmax_cell

  !--------------------------------------------------------------------!
  !                      retrieve neighbour list                       !
  !                                                                    !
  ! This routine returns the complete list of atom numbers of possible !
  ! neighbours.  It is not checked, if the returned atoms really are   !
  ! within the cut-off.                                                !
  !--------------------------------------------------------------------!

  subroutine lcl_nblist(iatom, nnb, nblist, stat)

    implicit none

    integer,                 intent(in)    :: iatom
    integer,                 intent(inout) :: nnb
    integer, dimension(nnb), intent(out)   :: nblist
    integer,       optional, intent(out)   :: stat

    integer, dimension(3) :: ic0, ic
    integer               :: inb, iat
    integer               :: iv, sgn

    if (.not. isInit) then
       write(0,*) "Error: module not initialized in `nblist'."
       stop
    end if

    if ((iatom < 1) .or. (iatom > nAtoms)) then
       write(0,*) "Error: invalid atom number in `nblist' : ", iatom
       stop
    end if

    if ((nnb < nmax_nblist) .and. (.not. present(stat))) then
       write(0,*) "Error: array nblist too small in `nblist'."
       stop
    end if

    if (present(stat)) stat = 0
    if (nnb == 0) return

    ! iatom is in cell ic0:
    ic0(1:3) = cellList(1:3,iatom)

    ! (1) other atoms in the same cell:

    inb = 0
    iat = cell(ic0(3),ic0(2),ic0(1))
    if (iat /= 0) then
       if (iat /= iatom) then
          inb = inb + 1
          if (inb > nnb) then
             if (present(stat)) stat = 1
             return
          end if
          nblist(inb) = iat
       end if
       do while(atomList(iat) /= 0)
          iat = atomList(iat)
          if (iat == iatom) cycle
          inb = inb + 1
          if (inb > nnb) then
             if (present(stat)) stat = 2
             return
          end if
          nblist(inb) = iat
       end do
    end if

    ! (2) atoms in neighbouring cells:

    do iv = 1, nCvecs
    do sgn = 1, -1, -2

       if (pbc) then
          ! cell coordinates obeying PBC
          ic(1) = mod(ic0(1) + sgn*Cvec(1,iv) + nCells(1) - 1, nCells(1)) + 1
          ic(2) = mod(ic0(2) + sgn*Cvec(2,iv) + nCells(2) - 1, nCells(2)) + 1
          ic(3) = mod(ic0(3) + sgn*Cvec(3,iv) + nCells(3) - 1, nCells(3)) + 1
       else
          ! isolated structures
          ic(1) = ic0(1) + sgn*Cvec(1,iv)
          ic(2) = ic0(2) + sgn*Cvec(2,iv)
          ic(3) = ic0(3) + sgn*Cvec(3,iv)
          if (any(ic > nCells) .or. any(ic < 1)) cycle
       end if

       ! all atoms in that cell:
       iat = cell(ic(3),ic(2),ic(1))
       if (iat /= 0) then
          if (.not. any(nblist(1:inb) == iat)) then
             inb = inb + 1
             if (inb > nnb) then
                if (present(stat)) stat = 3
                return
             end if
             nblist(inb) = iat
          end if
          do while(atomList(iat) /= 0)
             iat = atomList(iat)
             if (.not. any(nblist(1:inb) == iat)) then
                inb = inb + 1
                if (inb > nnb) then
                   if (present(stat)) stat = 4
                   return
                end if
                nblist(inb) = iat
             end if
          end do
       end if

    end do
    end do

    nnb = inb

  end subroutine lcl_nblist

  !--------------------------------------------------------------------!
  !      retrieve real PBC coordinates of the neighbouring atoms       !
  !--------------------------------------------------------------------!

  subroutine lcl_nbdist(iatom, nnb, nbcoo, nbdist, r_cut, itype, stat)

    implicit none

    integer,                            intent(in)    :: iatom
    integer,                            intent(inout) :: nnb
    double precision, dimension(3,nnb), intent(out)   :: nbcoo
    double precision, dimension(nnb),   intent(out)   :: nbdist
    double precision, optional,         intent(in)    :: r_cut
    integer,          optional,         intent(in)    :: itype
    integer,          optional,         intent(out)   :: stat

    integer,          dimension(nmax_nblist)          :: nblist
    integer                                           :: nnb_tot, inb, nnb2
    integer                                           :: iat, iT, sgn
    double precision                                  :: Rc, Rc2, dist2
    double precision, dimension(3)                    :: coo2, cart

    if (.not. isInit) then
       write(0,*) "Error: module not initialized in `nbdist'."
       stop
    end if

    if ((iatom < 1) .or. (iatom > nAtoms)) then
       write(0,*) "Error: invalid atom number in `nbdist' : ", iatom
       stop
    end if

    if ((nnb < nmax_nbdist) .and. (.not. present(stat))) then
       write(0,*) "Error: array nbcoo too small in `nbdist'. (1)"
       stop
    end if

    if (present(stat)) stat = 0
    if (present(r_cut)) then
       Rc = r_cut
    else
       Rc = r_max
    end if
    Rc2 = Rc*Rc
    nnb_tot = 0

    ! (1) get neighbour list:

    nnb2 = nmax_nblist
    call lcl_nblist(iatom, nnb2, nblist)

    ! (2) check distance do the periodic images of the central atom:

    nnb_tot = 0
    if ( (.not. present(itype)) .or. (atomType(iatom) == itype)) then
       do iT = 1, nTvecs
          if (TvecLen(iT) > Rc) exit
          nnb_tot = nnb_tot + 2
          if (nnb_tot > nnb) then
             if (present(stat)) then
                stat = 1
                return
             else
                write(0,*) "Error: array nbcoo too small in `nbdist'. (2)"
                stop
             end if
          end if ! nbcoo too small
          nbcoo(1,nnb_tot-1) = cooLatt(1,iatom) + dble(Tvec(1,iT))
          nbcoo(2,nnb_tot-1) = cooLatt(2,iatom) + dble(Tvec(2,iT))
          nbcoo(3,nnb_tot-1) = cooLatt(3,iatom) + dble(Tvec(3,iT))
          nbdist(nnb_tot-1)  = TvecLen(iT)
          nbcoo(1,nnb_tot)   = cooLatt(1,iatom) - dble(Tvec(1,iT))
          nbcoo(2,nnb_tot)   = cooLatt(2,iatom) - dble(Tvec(2,iT))
          nbcoo(3,nnb_tot)   = cooLatt(3,iatom) - dble(Tvec(3,iT))
          nbdist(nnb_tot)    = TvecLen(iT)
       end do
    end if ! correct atom type

    ! (3) calculate PBC distances for all atoms in the neighbourlist:

    do inb = 1, nnb2
       iat = nblist(inb)
       if ( present(itype) .and. (atomType(iat) /= itype)) cycle
       ! in home unit cell:
       cart(1:3) = cooLatt(1:3,iat) - cooLatt(1:3,iatom)
       cart(1:3) = matmul(latticeVec, cart(1:3))
       dist2 = sum(cart*cart)
       if (dist2 <= Rc2) then
          nnb_tot = nnb_tot + 1
          if (nnb_tot > nnb) then
             if (present(stat)) then
                stat = 2
                return
             else
                write(0,*) "Error: array nbcoo too small in `nbdist'. (3)"
                stop
             end if
          end if ! nbcoo too small
          nbcoo(1:3,nnb_tot) = cooLatt(1:3,iat)
          nbdist(nnb_tot)    = sqrt(dist2)
       end if ! within cut-off
       ! in periodic images:
       do iT = 1, nTvecs
       do sgn = 1, -1, -2
          coo2(1) = cooLatt(1,iat) + dble(sgn*Tvec(1,iT))
          coo2(2) = cooLatt(2,iat) + dble(sgn*Tvec(2,iT))
          coo2(3) = cooLatt(3,iat) + dble(sgn*Tvec(3,iT))
          cart(1:3) = coo2(1:3) - cooLatt(1:3,iatom)
          cart(1:3) = matmul(latticeVec, cart(1:3))
          dist2 = sum(cart*cart)
          if (dist2 <= Rc2) then
             nnb_tot = nnb_tot + 1
             if (nnb_tot > nnb) then
                if (present(stat)) then
                   stat = 3
                   return
                else
                   write(0,*) "Error: array nbcoo too small in `ndist'. (4)"
                   stop
                end if
             end if ! nbcoo too small
             nbcoo(1:3,nnb_tot) = coo2(1:3)
             nbdist(nnb_tot)    = sqrt(dist2)
          end if ! within cut-off
       end do
       end do
    end do

    ! (4) return number of neighbours found:

    nnb = nnb_tot

  end subroutine lcl_nbdist

  !--------------------------------------------------------------------!
  !    same as `lcl_nbdist', but cartesian coordinates are returned    !
  !--------------------------------------------------------------------!

  subroutine lcl_nbdist_cart(iatom, nnb, nbcoo, nbdist, r_cut, itype, &
                             nblist, nbtype, stat)

    implicit none

    integer,                                    intent(in)    :: iatom
    integer,                                    intent(inout) :: nnb
    double precision, dimension(3,nnb),         intent(out)   :: nbcoo
    double precision, dimension(nnb),           intent(out)   :: nbdist
    double precision,                 optional, intent(in)    :: r_cut
    integer,                          optional, intent(in)    :: itype
    integer,          dimension(nnb), optional, intent(out)   :: nblist
    integer,          dimension(nnb), optional, intent(out)   :: nbtype
    integer,                          optional, intent(out)   :: stat

    integer,          dimension(nmax_nblist)          :: nblist_loc
    integer                                           :: nnb_tot, inb, nnb2
    integer                                           :: iat, iT, sgn
    double precision                                  :: Rc, Rc2, dist2
    double precision, dimension(3)                    :: coo1, coo2, cart

    if (.not. isInit) then
       write(0,*) "Error: module not initialized in `nbdist'."
       stop
    end if

    if ((iatom < 1) .or. (iatom > nAtoms)) then
       write(0,*) "Error: invalid atom number in `nbdist' : ", iatom
       stop
    end if

    if ((nnb < nmax_nbdist) .and. (.not. present(stat))) then
       write(0,*) "Error: array nbcoo too small in `nbdist_cart'. (1)"
       stop
    end if

    if (present(stat)) stat = 0
    if (present(r_cut)) then
       Rc = r_cut
    else
       Rc = r_max
    end if
    Rc2 = Rc*Rc
    nnb_tot = 0

    ! (1) get neighbour list:

    nnb2 = nmax_nblist
    call lcl_nblist(iatom, nnb2, nblist_loc)

    ! (2) check distance do the periodic images of the central atom:

    nnb_tot = 0
    if ( (.not. present(itype)) .or. (atomType(iatom) == itype)) then
       do iT = 1, nTvecs
          if (TvecLen(iT) > Rc) exit
          nnb_tot = nnb_tot + 2
          if (nnb_tot > nnb) then
             if (present(stat)) then
                stat = 1
                return
             else
                write(0,*) "Error: array nbcoo too small in `nbdist_cart'. (2)"
                stop
             end if
          end if ! nbcoo too small
          if (present(nblist)) then
             nblist(nnb_tot-1) = iatom
             nblist(nnb_tot)   = iatom
          end if
          if (present(nbtype)) then
             nbtype(nnb_tot-1) = atomType(iatom)
             nbtype(nnb_tot)   = atomType(iatom)
          end if
          nbcoo(1,nnb_tot-1) = cooLatt(1,iatom) + dble(Tvec(1,iT))
          nbcoo(2,nnb_tot-1) = cooLatt(2,iatom) + dble(Tvec(2,iT))
          nbcoo(3,nnb_tot-1) = cooLatt(3,iatom) + dble(Tvec(3,iT))
          nbcoo(:,nnb_tot-1) = matmul(latticeVec, nbcoo(:,nnb_tot-1))
          nbdist(nnb_tot-1)  = TvecLen(iT)
          nbcoo(1,nnb_tot)   = cooLatt(1,iatom) - dble(Tvec(1,iT))
          nbcoo(2,nnb_tot)   = cooLatt(2,iatom) - dble(Tvec(2,iT))
          nbcoo(3,nnb_tot)   = cooLatt(3,iatom) - dble(Tvec(3,iT))
          nbcoo(:,nnb_tot)   = matmul(latticeVec, nbcoo(:,nnb_tot))
          nbdist(nnb_tot)    = TvecLen(iT)
       end do
    end if ! correct atom type

    ! (3) calculate PBC distances for all atoms in the neighbourlist:

    ! cartesian coordinates of central atom:
    coo1(1:3) = matmul(latticeVec, cooLatt(1:3,iatom))

    do inb = 1, nnb2
       iat = nblist_loc(inb)
       if ( present(itype) .and. (atomType(iat) /= itype)) cycle
       ! in home unit cell:
       coo2(1:3) = matmul(latticeVec, cooLatt(1:3,iat))
       cart(1:3) = coo2(1:3) - coo1(1:3)
       dist2 = sum(cart*cart)
       if (dist2 <= Rc2) then
          nnb_tot = nnb_tot + 1
          if (nnb_tot > nnb) then
             if (present(stat)) then
                stat = 2
                return
             else
                write(0,*) "Error: array nbcoo too small in `nbdist_cart'. (3)"
                stop
             end if
          end if ! nbcoo too small
          if(present(nblist)) nblist(nnb_tot) = iat
          if(present(nbtype)) nbtype(nnb_tot) = atomType(iat)
          nbcoo(1:3,nnb_tot) = coo2(1:3)
          nbdist(nnb_tot)    = sqrt(dist2)
       end if ! within cut-off
       ! in periodic images:
       do iT = 1, nTvecs
       do sgn = 1, -1, -2
          coo2(1) = cooLatt(1,iat) + dble(sgn*Tvec(1,iT))
          coo2(2) = cooLatt(2,iat) + dble(sgn*Tvec(2,iT))
          coo2(3) = cooLatt(3,iat) + dble(sgn*Tvec(3,iT))
          coo2(:) = matmul(latticeVec, coo2(1:3))
          cart(:) = coo2(1:3) - coo1(1:3)
          dist2 = sum(cart*cart)
          if (dist2 <= Rc2) then
             nnb_tot = nnb_tot + 1
             if (nnb_tot > nnb) then
                if (present(stat)) then
                   stat = 3
                   return
                else
                   write(0,*) "Error: array nbcoo too small in `nbdist_cart'. (4)"
                   stop
                end if
             end if ! nbcoo too small
             if(present(nblist)) nblist(nnb_tot) = iat
             if(present(nbtype)) nbtype(nnb_tot) = atomType(iat)
             nbcoo(1:3,nnb_tot) = coo2(1:3)
             nbdist(nnb_tot)    = sqrt(dist2)
          end if ! within cut-off
       end do
       end do
    end do

    ! (4) return number of neighbours found:

    nnb = nnb_tot

  end subroutine lcl_nbdist_cart


  !====================================================================!
  !                                                                    !
  !                         private procedures                         !
  !                                                                    !
  !====================================================================!



  !--------------------------------------------------------------------!
  !                       assign atoms to cells                        !
  !--------------------------------------------------------------------!

  subroutine cell_assign_atoms(nc, nAtoms, cooLatt, cell, atomList, &
                               cellList)

    implicit none

    integer,          dimension(3),                 intent(in)    :: nc
    integer,                                        intent(in)    :: nAtoms
    double precision, dimension(3,nAtoms),          intent(inout) :: cooLatt
    integer,          dimension(nc(3),nc(2),nc(1)), intent(out)   :: cell
    integer,          dimension(nAtoms),            intent(out)   :: atomList
    integer,          dimension(3,nAtoms),          intent(out)   :: cellList

    integer               :: iatom
    integer, dimension(3) :: ic

    cell(:,:,:)   = 0
    atomList(:)   = 0
    cellList(:,:) = 0

    do iatom = 1, nAtoms
       ! coordinates of isolated structures will be checked,
       ! but will not be wrapped
       call wrap_cooLatt(cooLatt(1:3,iatom))
       ic(1:3) = cell_get(cooLatt(1:3,iatom), nc)
       call cell_add(iatom, ic, nAtoms, nc, cell, atomList, cellList)
    end do

  end subroutine cell_assign_atoms

  !--------------------------------------------------------------------!
  !                determine the optimal cell division                 !
  !--------------------------------------------------------------------!

  subroutine cell_multiples(r_max, latticeVec, gridVec, nCells)

    implicit none

    double precision,                 intent(in)  :: r_max
    double precision, dimension(3,3), intent(in)  :: latticeVec
    double precision, dimension(3,3), intent(out) :: gridVec
    integer,          dimension(3),   intent(out) :: nCells

    double precision, dimension(3)   :: v
    double precision, dimension(3)   :: h
    double precision                 :: h_min
    double precision                 :: r

    ! compute the three heights:

    v  = vproduct(latticeVec(:,1),latticeVec(:,2))
    v  = v/sqrt(sum(v*v))
    h(3) = abs(sum(v*latticeVec(:,3)))

    v  = vproduct(latticeVec(:,1),latticeVec(:,3))
    v  = v/sqrt(sum(v*v))
    h(2) = abs(sum(v*latticeVec(:,2)))

    v  = vproduct(latticeVec(:,2),latticeVec(:,3))
    v  = v/sqrt(sum(v*v))
    h(1) = abs(sum(v*latticeVec(:,1)))

    ! shortest height:

    h_min = minval(h)

    ! compute multiples of h_min wrt. each height:

    nCells(1) = floor(h(1)/h_min)
    nCells(2) = floor(h(2)/h_min)
    nCells(3) = floor(h(3)/h_min)

    ! compute largest sphere that fits into such a cell:

    gridVec(:,1) = latticeVec(:,1)/dble(nCells(1))
    gridVec(:,2) = latticeVec(:,2)/dble(nCells(2))
    gridVec(:,3) = latticeVec(:,3)/dble(nCells(3))
    call inner_sphere(gridVec, r)

    ! reduce cell size if that sphere is unnecessary large:

    if (r > 0.5d0*r_max) then
       nCells(1) = floor(nCells(1)*(2.0d0*r/r_max))
       nCells(2) = floor(nCells(2)*(2.0d0*r/r_max))
       nCells(3) = floor(nCells(3)*(2.0d0*r/r_max))
       gridVec(:,1) = latticeVec(:,1)/dble(nCells(1))
       gridVec(:,2) = latticeVec(:,2)/dble(nCells(2))
       gridVec(:,3) = latticeVec(:,3)/dble(nCells(3))
    end if

  end subroutine cell_multiples


  !====================================================================!
  !                                                                    !
  !                 operations on the linked cell list                 !
  !                                                                    !
  !====================================================================!

  !--------------------------------------------------------------------!
  !             get cell for specific lattice coordinates              !
  !--------------------------------------------------------------------!

  function cell_get(coo, nc) result(ic)

    implicit none

    double precision, dimension(3), intent(in) :: coo
    integer,          dimension(3), intent(in) :: nc
    integer,          dimension(3)             :: ic

    ic(1) = nint(coo(1)*dble(nc(1)) + 0.5d0)
    ic(2) = nint(coo(2)*dble(nc(2)) + 0.5d0)
    ic(3) = nint(coo(3)*dble(nc(3)) + 0.5d0)

  end function cell_get

  !--------------------------------------------------------------------!
  !                          add atom to cell                          !
  !--------------------------------------------------------------------!

  subroutine cell_add(iatom, ic, nAtoms, nc, cell, atomList, cellList)

    implicit none

    integer,                               intent(in)    :: iatom
    integer, dimension(3),                 intent(in)    :: ic
    integer,                               intent(in)    :: nAtoms
    integer, dimension(3),                 intent(in)    :: nc
    integer, dimension(nc(3),nc(2),nc(1)), intent(inout) :: cell
    integer, dimension(nAtoms),            intent(inout) :: atomList
    integer, dimension(3,nAtoms),          intent(inout) :: cellList

    atomList(iatom) = cell(ic(3),ic(2),ic(1))
    cell(ic(3),ic(2),ic(1)) = iatom

    cellList(1:3,iatom) = ic(1:3)

  end subroutine cell_add

  !--------------------------------------------------------------------!
  !                       remove atom from cell                        !
  !--------------------------------------------------------------------!

!!$  subroutine cell_del(iatom, ic, nAtoms, nc, cell, atomList)
!!$
!!$    implicit none
!!$
!!$    integer,                               intent(in)    :: iatom
!!$    integer, dimension(3),                 intent(in)    :: ic
!!$    integer,                               intent(in)    :: nAtoms
!!$    integer, dimension(3),                 intent(in)    :: nc
!!$    integer, dimension(nc(3),nc(2),nc(1)), intent(inout) :: cell
!!$    integer, dimension(nAtoms),            intent(inout) :: atomList
!!$
!!$    integer :: i
!!$
!!$    if (cell(ic(3),ic(2),ic(1)) == iatom) then
!!$       cell(ic(3),ic(2),ic(1)) = atomList(iatom)
!!$    else
!!$       i = cell(ic(3),ic(2),ic(1))
!!$       do while(atomList(i) /= iatom)
!!$          if (i == 0) then
!!$             write(0,*) "Warning: atom not in cell in `cell_del'."
!!$             return
!!$          end if
!!$          i = atomList(i)
!!$       end do
!!$       atomList(i) = atomlist(iatom)
!!$    end if
!!$
!!$  end subroutine cell_del



  !====================================================================!
  !                                                                    !
  !                        auxilliary functions                        !
  !                                                                    !
  !====================================================================!



  !--------------------------------------------------------------------!
  !          wrap lattice coordinates to the interval [0,1[
  !--------------------------------------------------------------------!

  subroutine wrap_cooLatt(coo)

    implicit none

    double precision, dimension(3), intent(inout) :: coo

    if (.not. pbc) then
       if (any(coo < 0.0d0) .or. any(coo > 1.0d0)) then
          write(0,*) "Warning: atoms outside of bounding box in LC list!"
       end if
       return
    end if

    do while(coo(1) < 0.0d0)
       coo(1) = coo(1) + 1.0d0
    end do
    do while(coo(1) >= 1.0d0)
       coo(1) = coo(1) - 1.0d0
    end do

    do while(coo(2) < 0.0d0)
       coo(2) = coo(2) + 1.0d0
    end do
    do while(coo(2) >= 1.0d0)
       coo(2) = coo(2) - 1.0d0
    end do

    do while(coo(3) < 0.0d0)
       coo(3) = coo(3) + 1.0d0
    end do
    do while(coo(3) >= 1.0d0)
       coo(3) = coo(3) - 1.0d0
    end do

  end subroutine wrap_cooLatt

  !--------------------------------------------------------------------!
  !               volume of the cell spanned by a1,a2,a3               !
  !--------------------------------------------------------------------!

  function cell_volume(avec) result(V)

    implicit none

    double precision, dimension(3,3), intent(in) :: avec
    double precision                             :: V

    V = avec(1,1)*avec(2,2)*avec(3,3) &
      + avec(2,1)*avec(3,2)*avec(1,3) &
      + avec(3,1)*avec(1,2)*avec(2,3) &
      - avec(3,1)*avec(2,2)*avec(1,3) &
      - avec(1,1)*avec(3,2)*avec(2,3) &
      - avec(2,1)*avec(1,2)*avec(3,3)

    V = abs(V)

  end function cell_volume

  !--------------------------------------------------------------------!
  !          largest sphere within a cell spanned by a1,a2,a3          !
  !--------------------------------------------------------------------!

  subroutine inner_sphere(avec, r, c)

    implicit none

    double precision, dimension(3,3),           intent(in)  :: avec
    double precision,                           intent(out) :: r
    double precision, dimension(3),   optional, intent(out) :: c

    double precision, dimension(3) :: v
    double precision               :: h

    ! determine the shortest of the three heights:

    v = vproduct(avec(:,1),avec(:,2))
    v = v/sqrt(sum(v*v))
    h = abs(sum(v*avec(:,3)))
    r = h

    v = vproduct(avec(:,1),avec(:,3))
    v = v/sqrt(sum(v*v))
    h = abs(sum(v*avec(:,2)))
    r = min(r, h)

    v = vproduct(avec(:,2),avec(:,3))
    v = v/sqrt(sum(v*v))
    h = abs(sum(v*avec(:,1)))
    r = min(r, h)

    r = 0.5d0*r

    if (present(c)) c = 0.5d0*(avec(:,1) + avec(:,2) + avec(:,3))

  end subroutine inner_sphere

  !------------------------------------------------------------------!
  !                       vector/cross product                       !
  !------------------------------------------------------------------!

  function vproduct(a,b) result(c)

    implicit none

    double precision, dimension(3), intent(in) :: a, b
    double precision, dimension(3)             :: c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function vproduct

  !--------------------------------------------------------------------!
  !              periodic images within a cut-off radius               !
  !                                                                    !
  ! The subroutine returns the positive `half star' of translation     !
  ! vectors NOT including T = (0,0,0).                                 !
  !--------------------------------------------------------------------!

  subroutine translation_vectors(Rc, avec, nT, T, Tnorm, nc)

    implicit none

    !------------------------------------------------------------------!
    ! Rc         : cut-off radius                                      !
    ! avec(i,j)  : i-th component of lattice/cell vector j             !
    ! nT         : on entry: dimension of array T --> dimension(3,nT)  !
    !              on exit : number of T vectors found (can be > nT)   !
    !              --> the routine can be used just to calculate the   !
    !                  number of T vectors by setting nT = 0.          !
    ! T(i,j)     : (output) i-th component of the j-th T vector        !
    ! Tnorm(j)   : norm/length of the j-th T vector                    !
    !------------------------------------------------------------------!

    double precision,                          intent(in)    :: Rc
    double precision, dimension(3,3),          intent(in)    :: avec
    integer,                                   intent(inout) :: nT
    integer,          dimension(3,nT),         intent(out)   :: T
    double precision, dimension(nT), optional, intent(out)   :: Tnorm
    integer,          dimension(3),  optional, intent(in)    :: nc

    double precision,      parameter :: EPS = 1.0d-6

    double precision, dimension(3,8)              :: corner
    integer                                       :: ic, ic2, iT
    integer                                       :: i1, i2, i3
    integer                                       :: s2, s3
    logical                                       :: acc1, acc2, acc3
    double precision                              :: Rc2
    double precision, dimension(3)                :: v1, v2
    double precision                              :: vnorm1, vnorm2
    integer,          dimension(:),   allocatable :: idx
    integer,          dimension(:,:), allocatable :: T2
    double precision, dimension(:),   allocatable :: Tnorm2

    Rc2 = Rc*Rc
    iT  = 0

    ! 1/2 way to the corners of the cell spanned by avec:

    ic = 0
    do i3 = 0, 1
    do i2 = 0, 1
    do i1 = 0, 1
       ic = ic + 1
       corner(1:3,ic) = (dble(i1) - 0.5d0)*avec(1:3,1) &
                      + (dble(i2) - 0.5d0)*avec(1:3,2) &
                      + (dble(i3) - 0.5d0)*avec(1:3,3)
    end do
    end do
    end do

    ! loop over T vectors of increasing length:

    i3 = 0
    loop3 : do
       acc3 = .false.
       if (present(nc) .and. i3 >= nc(3)) exit loop3
       i2 = 0
       loop2 : do
          acc2 = .false.
          if (present(nc) .and. i2 >= nc(2)) exit loop2
          i1 = 0
          loop1 : do
             acc1 = .false.
             if (present(nc) .and. i1 >= nc(1)) exit loop1

             ! don't include T = (0,0,0)
             if (all((/i1,i2,i3/) == 0)) then
                acc2 = .true.
                i1 = i1 + 1
                cycle loop1
             end if

             sign3 : do s3 =  1, -1, -2
             sign2 : do s2 =  1, -1, -2

                v1(1:3) = dble(i1)*avec(1:3,1) &
                        + dble(s2*i2)*avec(1:3,2) &
                        + dble(s3*i3)*avec(1:3,3)

                ! loop over the eight corners of the cell:
                loopc  : do ic  = 1, 8
                loopc2 : do ic2 = 1, 8
                   v2(1:3) = v1(1:3) + corner(1:3,ic) - corner(1:3,ic2)
                   vnorm2 = sum(v2*v2)
                   if (vnorm2 - Rc2 < EPS) then
                      acc1 = .true.
                      iT = iT + 1
                      if (iT <= nT) then
                         T(1:3,iT) = (/ i1, s2*i2, s3*i3 /)
                         if (present(Tnorm)) then
                            vnorm1    = sum(v1*v1)
                            Tnorm(iT) = sqrt(vnorm1)
                         end if
                      end if
                      exit loopc
                   end if
                end do loopc2
                end do loopc

                if ((i1 == 0) .or. (i2 == 0)) exit sign2
             end do sign2
             if (((i1 == 0) .and. (i2 == 0)) .or. (i3 == 0)) exit sign3
             end do sign3

             if (.not. acc1) exit loop1
             acc2 = .true.
             i1 = i1 + 1
          end do loop1
          if (.not. acc2) exit loop2
          acc3 = .true.
          i2 = i2 + 1
       end do loop2
       if (.not. acc3) exit loop3
       i3 = i3 + 1
    end do loop3

    ! sort T vectors by length, but only if all of them were stored:

    if (present(Tnorm) .and. (iT <= nT)) then
       nT = iT
       allocate(idx(nT), T2(3,nT), Tnorm2(nT))
       call argsort(Tnorm(1:nT), idx)
       T2(:,1:nT)   = T(:,1:nT)
       Tnorm2(1:nT) = Tnorm(1:nT)
       do iT = 1, nT
          T(:,iT)   = T2(:,idx(iT))
          Tnorm(iT) = Tnorm2(idx(iT))
       end do
       deallocate(idx, T2, Tnorm2)
    else
       nT = iT
    end if

  end subroutine translation_vectors

  !--------------------------------------------------------------------!
  ! This procedure is not working properly, because it considers the   !
  ! outer sphere of the cell instead of the inner one.  Debug it when  !
  ! you feel bored.                                                    !
  !--------------------------------------------------------------------!

!!$  subroutine translation_vectors_old(Rc, avec, nT, T, Tnorm, nc)
!!$
!!$    !------------------------------------------------------------------!
!!$    ! Rc         : cut-off radius                                      !
!!$    ! avec(i,j)  : i-th component of lattice/cell vector j             !
!!$    ! nT         : on entry: dimension of array T --> dimension(3,nT)  !
!!$    !              on exit : number of T vectors found (can be > nT)   !
!!$    !              --> the routine can be used just to calculate the   !
!!$    !                  number of T vectors by setting nT = 0.          !
!!$    ! T(i,j)     : (output) i-th component of the j-th T vector        !
!!$    ! Tnorm(j)   : norm/length of the j-th T vector                    !
!!$    !------------------------------------------------------------------!
!!$
!!$    double precision,                          intent(in)    :: Rc
!!$    double precision, dimension(3,3),          intent(in)    :: avec
!!$    integer,                                   intent(inout) :: nT
!!$    integer,          dimension(3,nT),         intent(out)   :: T
!!$    double precision, dimension(nT), optional, intent(out)   :: Tnorm
!!$    integer,          dimension(3),  optional, intent(in)    :: nc
!!$
!!$    double precision,      parameter :: EPS = 1.0d-6
!!$
!!$    double precision, dimension(3,8)              :: corner
!!$    integer                                       :: ic, ic2
!!$    integer                                       :: i, i1, i2, i3
!!$    integer                                       :: s2, s3
!!$    double precision, dimension(3)                :: v, v2
!!$    double precision                              :: vnorm, v2norm
!!$    double precision, dimension(3)                :: a2
!!$    double precision                              :: f
!!$    integer,          dimension(3)                :: mult
!!$    integer                                       :: iT, nT_needed
!!$    integer,          dimension(:),   allocatable :: idx
!!$    integer,          dimension(:,:), allocatable :: T2
!!$    double precision, dimension(:),   allocatable :: Tnorm2
!!$
!!$    ! corners of the cell spanned by avec:
!!$
!!$    ic = 0
!!$    do i3 = 0, 1
!!$    do i2 = 0, 1
!!$    do i1 = 0, 1
!!$       ic = ic + 1
!!$       corner(1:3,ic) = dble(i1)*avec(1:3,1) &
!!$                      + dble(i2)*avec(1:3,2) &
!!$                      + dble(i3)*avec(1:3,3)
!!$    end do
!!$    end do
!!$    end do
!!$
!!$    ! find vector v(:) that connects the two corners that are
!!$    ! furthest away from each other:
!!$
!!$    v(:)  = corner(:,2) - corner(:,1)
!!$    vnorm = sum(v*v)
!!$    do ic = 1, 8
!!$    do ic2 = ic+1, 8
!!$       v2(:)  = corner(:,ic2) - corner(:,ic)
!!$       v2norm = sum(v2*v2)
!!$       if (vnorm < v2norm) then
!!$          vnorm = v2norm
!!$          v     = v2
!!$       end if
!!$    end do
!!$    end do
!!$    vnorm = sqrt(vnorm)
!!$
!!$    ! lengths^2 of the lattice vectors:
!!$
!!$    do i = 1, 3
!!$       a2(i) = sum(avec(:,i)*avec(:,i))
!!$    end do
!!$
!!$    ! Determine f_i for {f_i*a_i} i=1..3 span a supercell
!!$    ! that is large enough.  Round up to get full cells.
!!$    !
!!$    ! f_i = |v*a_i|/|a_i|^2 * (2Rc + |v|)/|v|
!!$
!!$    do i = 1, 3
!!$       f = abs(sum(v*avec(:,i)))
!!$       f = f/a2(i) * (2.0d0*Rc + vnorm)/vnorm
!!$       mult(i) = ceiling(0.5d0*(f - 1.0d0))
!!$    end do
!!$
!!$    ! In each lattice vector direction `i' we need -mult(i)...+mult(i)
!!$    ! cells.  Due to PBC, cells might be identical
!!$    ! --> make sure that mult(i) < nc(i)
!!$
!!$    if (present(nc)) then
!!$       do i = 1, 3
!!$          mult(i) = min(mult(i), nc(i)-1)
!!$       end do
!!$    end if
!!$
!!$    ! the number of T vectors (in the half star) are
!!$    ! nT_needed T vectors are needed
!!$
!!$    nT_needed = (2*mult(1) + 1)*(2*mult(2) + 1)*(2*mult(3) + 1)
!!$    nT_needed = (nT_needed - 1)/2
!!$    if (nT < nT_needed) then
!!$       ! array too small --> just return the number of T vectors
!!$       nT = nT_needed
!!$       return
!!$    end if
!!$
!!$    ! set up the T vectors:
!!$
!!$    iT = 0
!!$    do i1 = 0, mult(1)
!!$    do i2 = 0, mult(2)
!!$    do i3 = 0, mult(3)
!!$       if (all((/ i1, i2, i3/) == 0)) cycle
!!$       sgn2 : do s2 = 1, -1, -2
!!$       sgn3 : do s3 = 1, -1, -2
!!$          iT = iT + 1
!!$          T(:,iT) = (/ i1, s2*i2, s3*i3 /)
!!$          if (present(Tnorm)) then
!!$             v(1:3) =    i1*avec(:,1) &
!!$                    + s2*i2*avec(:,2) &
!!$                    + s3*i3*avec(:,3)
!!$             Tnorm(iT) = sqrt(sum(v*v))
!!$          end if
!!$          if (((i1 == 0) .and. (i2 == 0)) .or. (i3 == 0)) exit sgn3
!!$       end do sgn3
!!$       if ((i1 == 0) .or. (i2 == 0)) exit sgn2
!!$       end do sgn2
!!$    end do
!!$    end do
!!$    end do
!!$
!!$    ! sort T vectors by length, if their length was computed:
!!$
!!$    if (present(Tnorm)) then
!!$       nT = iT
!!$       allocate(idx(nT), T2(3,nT), Tnorm2(nT))
!!$       call argsort(Tnorm(1:nT), idx)
!!$       T2(:,1:nT)   = T(:,1:nT)
!!$       Tnorm2(1:nT) = Tnorm(1:nT)
!!$       do iT = 1, nT
!!$          T(:,iT)   = T2(:,idx(iT))
!!$          Tnorm(iT) = Tnorm2(idx(iT))
!!$       end do
!!$       deallocate(idx, T2, Tnorm2)
!!$    end if
!!$
!!$  end subroutine translation_vectors_old

end module lclist
