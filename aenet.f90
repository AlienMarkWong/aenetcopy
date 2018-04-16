!-----------------------------------------------------------------------
! aenet.f90 -- public interface to the aenet routines (aenetLib)
!
! This module also provides the C bindings for libaenet.so/.a .
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
! 2014-08-31 Nongnuch Artrith, Alexander Urban
!-----------------------------------------------------------------------
module aenet

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_char, c_bool, &
                                         c_ptr, c_f_pointer, c_null_char

  use aeio,        only: TYPELEN, PATHLEN

  use feedforward, only: ff_eval, ff_deriv

  use geometry,    only: geo_type_conv

  use io,          only: io_cstring_len,      &
                         io_cstring2f

  use lclist,      only: lcl_nmax_nbdist

  use potential,   only: NNPot,               &
                         load_NNPot,          &
                         del_NNPot,           &
                         pot_init,            &
                         pot_final,           &
                         pot_print_info,      &
                         pot_get_range

  use sfsetup,     only: stp_init,            &
                         stp_final,           &
                         stp_nsf_max,         &
                         stp_eval,            &
                         sfval, sfderiv_i, sfderiv_j

  use timing,      only: tng_timing, tng_timing2

  implicit none
  private
  save

  public :: aenet_init,                     &
            aenet_final,                    &
            aenet_all_loaded,               &
            aenet_atomic_energy,            &
            aenet_atomic_energy_and_forces, &
            aenet_convert_atom_types,       &
            aenet_free_atom_energy,         &
            aenet_load_potential,           &
            aenet_print_info

  !---------------------------- constants -----------------------------!

  ! return status
  integer(kind=c_int), bind(C, name='AENET_OK'),         public :: AENET_OK = 0_c_int
  integer(kind=c_int), bind(C, name='AENET_ERR_INIT'),   public :: AENET_ERR_INIT = 1_c_int
  integer(kind=c_int), bind(C, name='AENET_ERR_MALLOC'), public :: AENET_ERR_MALLOC = 2_c_int
  integer(kind=c_int), bind(C, name='AENET_ERR_IO'),     public :: AENET_ERR_IO = 3_c_int
  integer(kind=c_int), bind(C, name='AENET_ERR_TYPE'),   public :: AENET_ERR_TYPE = 4_c_int

  integer(kind=c_int), bind(C, name='AENET_TYPELEN'),    public :: AENET_TYPELEN = TYPELEN
  integer(kind=c_int), bind(C, name='AENET_PATHLEN'),    public :: AENET_PATHLEN = PATHLEN

  !---------------------------- variables -----------------------------!

  integer(kind=c_int), bind(C, name='aenet_nsf_max'), public :: aenet_nsf_max
  integer(kind=c_int), bind(C, name='aenet_nnb_max'), public :: aenet_nnb_max
  real(kind=c_double), bind(C, name='aenet_Rc_min'),  public :: aenet_Rc_min
  real(kind=c_double), bind(C, name='aenet_Rc_max'),  public :: aenet_Rc_max

  !----------------------------- private ------------------------------!

  logical, private :: aenet_is_init = .false.
  logical, private :: aenet_is_loaded = .false.

  integer,                                           private :: aenet_ntypes
  character(len=TYPELEN), dimension(:), allocatable, private :: aenet_atom_types
  type(NNPot),            dimension(:), allocatable, private :: aenet_pot
  double precision,       dimension(:), allocatable, private :: aenet_dE_dG

contains

  !--------------------------------------------------------------------!
  !                  Initialization and Finalization                   !
  !--------------------------------------------------------------------!

  subroutine aenet_init(atom_types, stat)

    implicit none

    character(len=*), dimension(:), intent(in)  :: atom_types
    integer,                        intent(out) :: stat

    integer :: ok

    stat = AENET_OK
    if (aenet_is_init) then
       stat = AENET_ERR_INIT
       return
    end if

    aenet_ntypes = size(atom_types)
    allocate(aenet_pot(aenet_ntypes),        &
             aenet_atom_types(aenet_ntypes), &
             stat=ok)
    if (ok /= 0) then
       aenet_ntypes = 0
       stat = AENET_ERR_MALLOC
       return
    end if
    aenet_atom_types = atom_types
    aenet_is_init = .true.

  end subroutine aenet_init

  subroutine aenet_init_C(ntypes, atom_types, stat) bind(C, name='aenet_init')

    implicit none

    integer(kind=c_int), value,             intent(in)  :: nTypes
    type(c_ptr), dimension(ntypes), target, intent(in)  :: atom_types
    integer(kind=c_int),                    intent(out) :: stat

    character,                dimension(:), pointer :: fstringp
    character(len=TYPELEN+1), dimension(ntypes)     :: fatom_types

    integer :: i, slen

    do i = 1, ntypes
       call c_f_pointer(atom_types(i), fstringp, [TYPELEN+1])
       fatom_types(i) = transfer(fstringp(1:TYPELEN+1), fatom_types(i))
       slen = index(fatom_types(i), c_null_char) - 1
       fatom_types(i)(slen+1:TYPELEN+1) = ' '
    end do

    call aenet_init(fatom_types, stat)

  end subroutine aenet_init_C

  !--------------------------------------------------------------------!

  subroutine aenet_final(stat) bind(C)

    implicit none

    integer(kind=c_int), intent(out) :: stat

    integer :: itype
    integer :: ok

    stat = AENET_OK
    if (aenet_is_init) then
       call stp_final(aenet_ntypes, aenet_pot(1:aenet_ntypes)%stp)
       do itype = 1, aenet_ntypes
          call del_NNPot(aenet_pot(itype))
       end do
       deallocate(aenet_pot, aenet_atom_types, stat=ok)
       if (ok /= 0) then
          stat = AENET_ERR_MALLOC
          return
       end if
       aenet_ntypes = 0
       aenet_is_init = .false.
       aenet_is_loaded = .false.
    end if

    if (allocated(aenet_dE_dG)) deallocate(aenet_dE_dG)

  end subroutine aenet_final

  !--------------------------------------------------------------------!
  !                               output                               !
  !--------------------------------------------------------------------!

  subroutine aenet_print_info() bind(C)

    implicit none

    integer :: ipot

    if (.not. aenet_is_init) then
       write(*,*) "Nothing to print. AenetLib is not initialized."
    else
       do ipot = 1, aenet_ntypes
          if (aenet_pot(ipot)%init) then
             call pot_print_info(aenet_pot(ipot))
          end if
       end do
    end if

  end subroutine aenet_print_info

  !--------------------------------------------------------------------!
  !                        Load ANN potentials                         !
  !--------------------------------------------------------------------!

  subroutine aenet_load_potential(type_id, filename, stat)

    implicit none

    integer,          intent(in)  :: type_id
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: stat

    integer :: ok
    logical :: fexists

    stat = AENET_OK
    if (.not. aenet_is_init) then
       stat = AENET_ERR_INIT
       return
    end if

    if ((type_id <= 0) .or. (type_id > aenet_ntypes)) then
       stat = AENET_ERR_TYPE
       return
    end if

    inquire(file=trim(filename), exist=fexists)
    if (.not. fexists) then
       stat = AENET_ERR_IO
       return
    end if

    aenet_pot(type_id) = load_NNPot(aenet_atom_types, filename)

    ! when all potentials are loaded, determine array sizes
    if (aenet_all_loaded()) then
       call pot_get_range(aenet_ntypes, aenet_pot, aenet_Rc_min, aenet_Rc_max)
       if (stat /= 0) return
       aenet_nnb_max = lcl_nmax_nbdist(aenet_Rc_min, aenet_Rc_max)
       call stp_init(aenet_ntypes, aenet_pot(1:aenet_ntypes)%stp, aenet_nnb_max)
       aenet_nsf_max = stp_nsf_max()
       allocate(aenet_dE_dG(aenet_nsf_max), stat=ok)
       if (ok /= 0) then
          stat = AENET_ERR_MALLOC
          return
       end if
       aenet_is_loaded = .true.
    end if

  end subroutine aenet_load_potential

  subroutine aenet_load_potential_C(type_id, filename, stat) &
       bind(C, name='aenet_load_potential')

    implicit none

    integer(kind=c_int), value,           intent(in)  :: type_id
    character(kind=c_char), dimension(*), intent(in)  :: filename
    integer(kind=c_int),                  intent(out) :: stat

    integer :: slen
    character(len=1024) :: ffilename

    slen = io_cstring_len(filename)
    if (slen > 1024) then
       stat = aenet_ERR_MALLOC
       return
    end if
    ffilename = io_cstring2f(filename, slen)

    call aenet_load_potential(type_id, trim(ffilename), stat)

  end subroutine aenet_load_potential_C

  function aenet_all_loaded() result(all_loaded) bind(C)

    implicit none

    logical(kind=c_bool) :: all_loaded
    integer :: ipot

    all_loaded = .true.
    do ipot = 1, aenet_ntypes
       if (.not. aenet_pot(ipot)%init) then
          all_loaded = .false.
          return
       end if
    end do

  end function aenet_all_loaded

  !--------------------------------------------------------------------!
  !                            information                             !
  !--------------------------------------------------------------------!

  function aenet_free_atom_energy(type_id) result(E_atom) bind(C)

    implicit none

    integer(kind=c_int), intent(in) :: type_id
    real(kind=c_double)             :: E_atom

    E_atom = aenet_pot(type_id)%E_atom

  end function aenet_free_atom_energy

  !--------------------------------------------------------------------!
  !                             Evaluation                             !
  !                                                                    !
  ! Attention: all routines require synchronized atom type IDs, i.e.,  !
  !            the IDs passed to the evaluation routines must be       !
  !            compatible with the ANN potential type IDs.             !
  !                                                                    !
  ! Notes:     * Coordinates are Cartesian.                            !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine aenet_atomic_energy(coo_i, type_i, n_j, coo_j, type_j, &
                                 E_i, stat) bind(C)

    implicit none

    real(kind=c_double), dimension(3),     intent(in)  :: coo_i
    integer(kind=c_int), value,            intent(in)  :: type_i
    integer(kind=c_int), value,            intent(in)  :: n_j
    real(kind=c_double), dimension(3,n_j), intent(in)  :: coo_j
    integer(kind=c_int), dimension(n_j),   intent(in)  :: type_j
    real(kind=c_double),                   intent(out) :: E_i
    integer(kind=c_int),                   intent(out) :: stat

    double precision, dimension(1) :: E_i_arr
    integer :: nsf

    stat = AENET_OK
    if (.not. (aenet_is_init .and. aenet_is_loaded)) then
       stat = aenet_ERR_INIT
       return
    end if

    nsf = aenet_pot(type_i)%stp%nsf
    call stp_eval(type_i, coo_i, n_j, coo_j, type_j, &
                  aenet_pot(type_i)%stp, scaled=.true.)
    call ff_eval(aenet_pot(type_i)%net, nsf, sfval, 1, E_i_arr)
    E_i = aenet_pot(type_i)%E_scale*E_i_arr(1) + aenet_pot(type_i)%E_shift

    E_i = E_i + aenet_pot(type_i)%E_atom

  end subroutine aenet_atomic_energy

  subroutine aenet_atomic_energy_and_forces( &
       coo_i, type_i, index_i, n_j, coo_j, type_j, index_j, natoms, &
       E_i, F, stat) bind(C)

    implicit none

    real(kind=c_double), dimension(3),        intent(in)    :: coo_i
    integer(kind=c_int), value,               intent(in)    :: type_i
    integer(kind=c_int), value,               intent(in)    :: index_i
    integer(kind=c_int), value,               intent(in)    :: n_j
    real(kind=c_double), dimension(3,n_j),    intent(in)    :: coo_j
    integer(kind=c_int), dimension(n_j),      intent(in)    :: type_j
    integer(kind=c_int), dimension(n_j),      intent(in)    :: index_j
    integer(kind=c_int), value,               intent(in)    :: natoms
    real(kind=c_double),                      intent(out)   :: E_i
    real(kind=c_double), dimension(3,natoms), intent(inout) :: F
    integer(kind=c_int),                      intent(out)   :: stat

    double precision, dimension(1)              :: E_i_arr
    integer                                     :: nsf, j

    stat = AENET_OK
    if (.not. (aenet_is_init .and. aenet_is_loaded)) then
       stat = aenet_ERR_INIT
       return
    end if

    nsf = aenet_pot(type_i)%stp%nsf
    call stp_eval(type_i, coo_i, n_j, coo_j, type_j, &
                  aenet_pot(type_i)%stp, scaled=.true., deriv=.true.)

    call ff_eval(aenet_pot(type_i)%net, nsf, sfval, 1, E_i_arr)
    call ff_deriv(aenet_pot(type_i)%net, nsf, 1, aenet_dE_dG(1:nsf))

    E_i = aenet_pot(type_i)%E_scale*E_i_arr(1) + aenet_pot(type_i)%E_shift
    E_i = E_i + aenet_pot(type_i)%E_atom

    F(1:3, index_i) = F(1:3, index_i) - aenet_pot(type_i)%E_scale &
                    * matmul(sfderiv_i(1:3,1:nsf), aenet_dE_dG(1:nsf))

    do j = 1, n_j
       F(1:3, index_j(j)) = F(1:3, index_j(j)) - aenet_pot(type_i)%E_scale &
                          * matmul(sfderiv_j(1:3,1:nsf,j), aenet_dE_dG(1:nsf))
    end do

  end subroutine aenet_atomic_energy_and_forces

  !--------------------------------------------------------------------!
  !                       convert atom type IDs                        !
  !--------------------------------------------------------------------!

  subroutine aenet_convert_atom_types(&
       atom_types_in, type_id_in, type_id_out, stat)

    implicit none

    character(len=*), dimension(:), intent(in)  :: atom_types_in
    integer,          dimension(:), intent(in)  :: type_id_in
    integer,          dimension(:), intent(out) :: type_id_out
    integer,                        intent(out) :: stat

    integer :: iat,it,  nTypes_in, natoms_in

    stat = AENET_OK
    if (.not. aenet_is_init) then
       stat = aenet_ERR_INIT
       return
    end if

    ntypes_in = size(atom_types_in)
    natoms_in = size(type_id_in)

    do iat = 1, natoms_in
       it = geo_type_conv(type_id_in(iat), ntypes_in, atom_types_in, &
                          aenet_ntypes, aenet_atom_types)
       if (it == 0) then
          stat = aenet_ERR_TYPE
          return
       end if
       type_id_out(iat) = it
    end do

  end subroutine aenet_convert_atom_types

  subroutine aenet_convert_atom_types_C(ntypes_in, atom_types_in, &
                         natoms_in, type_id_in, type_id_out, stat &
                         ) bind(C, name="aenet_convert_atom_types")

    implicit none

    integer(kind=c_int), value,                        intent(in)  :: ntypes_in
    type(c_ptr),         dimension(ntypes_in), target, intent(in)  :: atom_types_in
    integer(kind=c_int), value,                        intent(in)  :: natoms_in
    integer(kind=c_int), dimension(natoms_in),         intent(in)  :: type_id_in
    integer(kind=c_int), dimension(natoms_in),         intent(out) :: type_id_out
    integer(kind=c_int),                               intent(out) :: stat

    character,                dimension(:), pointer :: fstringp
    character(len=TYPELEN+1), dimension(ntypes_in)  :: fatom_types

    integer :: i, slen

    do i = 1, ntypes_in
       call c_f_pointer(atom_types_in(i), fstringp, [TYPELEN+1])
       fatom_types(i) = transfer(fstringp(1:TYPELEN+1), fatom_types(i))
       slen = index(fatom_types(i), c_null_char) - 1
       fatom_types(i)(slen+1:TYPELEN+1) = ' '
    end do

    call aenet_convert_atom_types(&
         fatom_types, type_id_in, type_id_out, stat)

  end subroutine aenet_convert_atom_types_C

  !--------------------------------------------------------------------!
  !                         internal routines                          !
  !--------------------------------------------------------------------!

end module aenet
