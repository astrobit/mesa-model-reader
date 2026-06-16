! ==============================================================================
! MODULE: mesa_iso_mass
!
! Provides data structures and routines for accumulating the integrated mass
! of every isotope (and parent element) present in the MESA reaction network
! mesa_3335.net.
!
! Typical workflow
! ----------------
!   1. call init_iso_mass()          ! populate tables; zero all masses
!   2. call add_iso_mass(id, dm)     ! add a shell mass contribution dm [g]
!      ... (once per isotope per stellar zone, inside your integration loop)
!   3. Read isotope_table(:)%mass    ! integrated isotope masses [g]
!      Read element_table(:)%mass    ! integrated element masses [g]
!   4. call reset_iso_mass()         ! zero masses without re-building tables
!
! Data model
! ----------
!   isotope_entry  – one isotope: ID string, Z, A, accumulated mass [g]
!   element_entry  – one element: symbol, Z, accumulated mass [g]
!
! Isotope ID strings are constructed as  Cc<A>  where Cc is the title-cased
! chemical symbol and A is the mass number, e.g. "H1", "He4", "Fe56".
! These match the naming convention used by mesa_mod_summary.f90.
!
! Network coverage (from mesa_3335.net)
! --------------------------------------
! The network spans 89 species groups from the free neutron through Astatine
! (At).  Each element contributes all isotopes from its minimum to maximum
! mass number as listed in the .net file.  The total isotope count is 3335
! (matching the network name).
!
! Performance note
! ----------------
! add_iso_mass uses a binary search on the sorted isotope_table so that
! the per-zone overhead is O(log N) rather than O(N).
! ==============================================================================
module mesa_iso_mass
  use iso_fortran_env, only: real64
  implicit none
  private


  type, public :: element_type
    character(len=3) :: symbol = ''   ! ISO symbol, e.g. "Fe"
    integer          :: Z      = 0    ! atomic number
  end type element_type

    ! --------------------------------------------------------------------------
    ! Element specification table 
    ! Columns: ISO element symbol 
    !
    ! --------------------------------------------------------------------------
    integer, parameter :: N_SPEC = 119   ! number of element 

    type(element_type), parameter :: elem_table(N_SPEC) = [ &
      element_type('n  ',0), &   !  Z= 0  (free neutron,  A=1)
      element_type('H  ',1), &   !  Z= 1   
      element_type('He ',2), &   !  Z= 2   
      element_type('Li ',3), &   !  Z= 3   
      element_type('Be ',4), &   !  Z= 4   
      element_type('B  ',5), &   !  Z= 5   
      element_type('C  ',6), &   !  Z= 6   
      element_type('N  ',7), &   !  Z= 7   
      element_type('O  ',8), &   !  Z= 8   
      element_type('F  ',9), &   !  Z= 9   
      element_type('Ne ',10), &   !  Z=10   
      element_type('Na ',11), &   !  Z=11   
      element_type('Mg ',12), &   !  Z=12   
      element_type('Al ',13), &   !  Z=13   
      element_type('Si ',14), &   !  Z=14   
      element_type('P  ',15), &   !  Z=15   
      element_type('S  ',16), &   !  Z=16   
      element_type('Cl ',17), &   !  Z=17   
      element_type('Ar ',18), &   !  Z=18   
      element_type('K  ',19), &   !  Z=19   
      element_type('Ca ',20), &   !  Z=20   
      element_type('Sc ',21), &   !  Z=21   
      element_type('Ti ',22), &   !  Z=22   
      element_type('V  ',23), &   !  Z=23   
      element_type('Cr ',24), &   !  Z=24   
      element_type('Mn ',25), &   !  Z=25   
      element_type('Fe ',26), &   !  Z=26   
      element_type('Co ',27), &   !  Z=27   
      element_type('Ni ',28), &   !  Z=28   
      element_type('Cu ',29), &   !  Z=29   
      element_type('Zn ',30), &   !  Z=30   
      element_type('Ga ',31), &   !  Z=31   
      element_type('Ge ',32), &   !  Z=32   
      element_type('As ',33), &   !  Z=33   
      element_type('Se ',34), &   !  Z=34   
      element_type('Br ',35), &   !  Z=35   
      element_type('Kr ',36), &   !  Z=36
      element_type('Rb ',37), &   !  Z=37
      element_type('Sr ',38), &   !  Z=38
      element_type('Y  ',39), &   !  Z=39
      element_type('Zr ',40), &   !  Z=40
      element_type('Nb ',41), &   !  Z=41
      element_type('Mo ',42), &   !  Z=42
      element_type('Tc ',43), &   !  Z=43
      element_type('Ru ',44), &   !  Z=44
      element_type('Rh ',45), &   !  Z=45
      element_type('Pd ',46), &   !  Z=46
      element_type('Ag ',47), &   !  Z=47
      element_type('Cd ',48), &   !  Z=48
      element_type('In ',49), &   !  Z=49
      element_type('Sn ',50), &   !  Z=50
      element_type('Sb ',51), &   !  Z=51
      element_type('Te ',52), &   !  Z=52
      element_type('I  ',53), &   !  Z=53
      element_type('Xe ',54), &   !  Z=54
      element_type('Cs ',55), &   !  Z=55
      element_type('Ba ',56), &   !  Z=56
      element_type('La ',57), &   !  Z=57
      element_type('Ce ',58), &   !  Z=58
      element_type('Pr ',59), &   !  Z=59
      element_type('Nd ',60), &   !  Z=60
      element_type('Pm ',61), &   !  Z=61
      element_type('Sm ',62), &   !  Z=62
      element_type('Eu ',63), &   !  Z=63
      element_type('Gd ',64), &   !  Z=64
      element_type('Tb ',65), &   !  Z=65
      element_type('Dy ',66), &   !  Z=66
      element_type('Ho ',67), &   !  Z=67
      element_type('Er ',68), &   !  Z=68
      element_type('Tm ',69), &   !  Z=69
      element_type('Yb ',70), &   !  Z=70
      element_type('Lu ',71), &   !  Z=71
      element_type('Hf ',72), &   !  Z=72
      element_type('Ta ',73), &   !  Z=73
      element_type('W  ',74), &   !  Z=74
      element_type('Re ',75), &   !  Z=75
      element_type('Os ',76), &   !  Z=76
      element_type('Ir ',77), &   !  Z=77
      element_type('Pt ',78), &   !  Z=78
      element_type('Au ',79), &   !  Z=79
      element_type('Hg ',80), &   !  Z=80
      element_type('Tl ',81), &   !  Z=81
      element_type('Pb ',82), &   !  Z=82
      element_type('Bi ',83), &   !  Z=83
      element_type('Po ',84), &   !  Z=84
      element_type('At ',85), &   !  Z=85
      element_type('Rn ',86), &   !  Z=86
      element_type('Fr ',87), &   !  Z=87
      element_type('Ra ',88), &   !  Z=88
      element_type('Ac ',89), &   !  Z=89
      element_type('Th ',90), &   !  Z=90
      element_type('Pa ',91), &   !  Z=91
      element_type('U  ',92), &   !  Z=92      
      element_type('Np ',93), &   !  Z=93
      element_type('Pu ',94), &   !  Z=94
      element_type('Am ',95), &   !  Z=95
      element_type('Cm ',96), &   !  Z=96
      element_type('Bk ',97), &   !  Z=97
      element_type('Cf ',98), &   !  Z=98
      element_type('Es ',99), &   !  Z=99
      element_type('Fm ',100), &   !  Z=100
      element_type('Md ',101), &   !  Z=101
      element_type('No ',102), &   !  Z=102
      element_type('Lr ',103), &   !  Z=103
      element_type('Rf ',104), &   !  Z=104
      element_type('Db ',105), &   !  Z=105
      element_type('Sg ',106), &   !  Z=106
      element_type('Bh ',107), &   !  Z=107
      element_type('Hs ',108), &   !  Z=108
      element_type('Mt ',109), &   !  Z=109
      element_type('Ds ',110), &   !  Z=110
      element_type('Rg ',111), &   !  Z=111
      element_type('Cn ',112), &   !  Z=112
      element_type('Nh ',113), &   !  Z=113
      element_type('Fl ',114), &   !  Z=114
      element_type('Mc ',115), &   !  Z=115
      element_type('Lc ',116), &   !  Z=116
      element_type('Ts ',117), &   !  Z=117
      element_type('Og ',118) &   !  Z=118
    ]

  ! ---------------------------------------------------------------------------
  ! Public interface
  ! ---------------------------------------------------------------------------
  public :: isotope_entry
  public :: element_entry
  public :: isotope_table
  public :: element_table
  public :: n_isotopes
  public :: n_elements
  public :: init_iso_mass
  public :: reset_iso_mass
  public :: add_iso_mass
  public :: get_iso_mass
  public :: get_iso_mass_by_idx
  public :: find_element_by_Z
  public :: get_ele_mass
  public :: get_ele_mass_by_Z
  public :: title_case

  ! ---------------------------------------------------------------------------
  ! isotope_entry: one row of the isotope table
  !
  !   id      – title-cased symbol + mass number, e.g. "Fe56"
  !   symbol  – lower-case element symbol as it appears in MESA, e.g. "fe"
  !   Z       – atomic number (proton count)
  !   A       – mass number  (proton + neutron count)
  !   mass    – double-precision accumulated mass [g]; zero until add_iso_mass
  !             is called
  ! ---------------------------------------------------------------------------
  type, public :: isotope_entry
    character(len=6) :: id     = ''   ! e.g. "Fe56"
    type(element_type) :: element = element_type("",-1)
    integer          :: A      = 0    ! mass number
    real(real64)     :: mass   = 0.0d0
  end type isotope_entry

  ! ---------------------------------------------------------------------------
  ! element_entry: one row of the element table
  !
  !   symbol  – title-cased element symbol, e.g. "Fe"
  !   Z       – atomic number
  !   mass    – sum of all isotope masses for this element [g]
  ! ---------------------------------------------------------------------------
  type, public :: element_entry
    type(element_type) :: element = element_type("",-1)
    real(real64)     :: mass   = 0.0d0
  end type element_entry

  ! ---------------------------------------------------------------------------
  ! Module-level tables (allocated and populated by init_iso_mass)
  ! ---------------------------------------------------------------------------
  type(isotope_entry), allocatable :: isotope_table(:)
  type(element_entry), allocatable :: element_table(:)
  integer :: n_isotopes = 0
  integer :: n_elements = 0

  ! Private flag: tables have been built at least once
  logical, private :: tables_built = .false.


contains


  ! ---------------------------------------------------------------------------
  ! find_element_by_Z
  !
  ! Returns the index into element_table whose Z matches the argument,
  ! or 0 if not found.
  ! ---------------------------------------------------------------------------
  function find_element_by_Z(Z) result(elem)
    integer, intent(in) :: Z
    type(element_type) :: elem


    if (Z < N_SPEC) then
      elem = elem_table(Z + 1); ! Z+1 because the table include neutrons
    else
      write(0,'(3a)') 'mesa_iso_mass WARNING: find_element_symbol_by_Z: Z = ', Z, ' not found.'
    end if
    return
  end function find_element_by_Z


  ! ===========================================================================
  ! init_iso_mass
  !
  ! Builds isotope_table and element_table from the hardcoded mesa_3335.net
  ! isotope ranges and zeroes all accumulated masses.  Safe to call multiple
  ! times; calling it again re-zeroes the masses and rebuilds the tables.
  !
  ! The element specification block below mirrors mesa_3335.net exactly.
  ! Each row is:  (symbol, Z, A_min, A_max)
  ! A singleton isotope (e.g. "neut") uses A_min == A_max.
  ! ===========================================================================
  subroutine init_iso_mass()

    integer, parameter :: NISO = 86
    
    type :: isotope_specification
      integer :: Z
      integer :: Amin
      integer :: Amax
    end type isotope_specification
    
    type(isotope_specification), parameter :: isotope_specs(NISO) = [ &
      isotope_specification(0, 1, 1), & !n
      isotope_specification(1, 1, 3), & !H
      isotope_specification(2, 3, 4), & !He
      isotope_specification(3, 6, 9), &  ! lithium
      isotope_specification(4, 7, 10), & ! berylium
      isotope_specification(5, 8, 14), & ! boron
      isotope_specification(6, 9, 18), & ! carbon
      isotope_specification(7, 11, 21), & ! nitrogen
      isotope_specification(8, 13, 22), & ! oxygen
      isotope_specification(9, 16, 26), & ! fluorine
      isotope_specification(10, 16, 31), & ! neon
      isotope_specification(11, 19, 34), & ! sodium
      isotope_specification(12, 20, 37), & ! magnesium
      isotope_specification(13, 22, 40), & ! aluminum
      isotope_specification(14, 22, 43), & ! silicon
      isotope_specification(15, 24, 46), & !phosphorus
      isotope_specification(16, 26, 49), & !sulfer
      isotope_specification(17, 28, 51), & !clorine
      isotope_specification(18, 30, 54), & !argon
      isotope_specification(19, 32, 56), & !potasium
      isotope_specification(20, 34, 59), & !calcium
      isotope_specification(21, 36, 64), & !scandium
      isotope_specification(22, 38, 67), & !titanium
      isotope_specification(23, 40, 72), & !vandium
      isotope_specification(24, 42, 75), & !chromium
      isotope_specification(25, 44, 76), & !manganese
      isotope_specification(26, 46, 78), & !iron
      isotope_specification(27, 48, 80), & !cobolt
      isotope_specification(28, 50, 83), & !nickel
      isotope_specification(29, 52, 86), & !copper
      isotope_specification(30, 54, 89), & !zinc
      isotope_specification(31, 56, 92), & !gallium
      isotope_specification(32, 58, 95), & !germanium
      isotope_specification(33, 60, 100), & !arsenic
      isotope_specification(34, 63, 103), & !selenium
      isotope_specification(35, 65, 105), & !bromine
      isotope_specification(36, 68, 108), & !krypton
      isotope_specification(37, 73, 111), & !rubidium
      isotope_specification(38, 72, 114), & !strontium
      isotope_specification(39, 75, 119), & !yttrium
      isotope_specification(40, 76, 122), & !zirconium
      isotope_specification(41, 80, 124), & !niobium
      isotope_specification(42, 80, 125), & !molybdenum
      isotope_specification(43, 85, 126), & !technetium
      isotope_specification(44, 86, 128), & !ruthenium
      isotope_specification(45, 89, 130), & !rhodium
      isotope_specification(46, 88, 132), & !palladium
      isotope_specification(47, 93, 134), & !silver
      isotope_specification(48, 95, 134), & !cadmium
      isotope_specification(49, 97, 149), & !indium
      isotope_specification(50, 99, 152), & !tin
      isotope_specification(51, 105, 151), & !antimony
      isotope_specification(52, 104, 154), & !tellurium
      isotope_specification(53, 108, 158), & !iodine
      isotope_specification(54, 108, 161), & !xenon
      isotope_specification(55, 118, 165), & !cesium
      isotope_specification(56, 119, 168), & !barium
      isotope_specification(57, 122, 172), & !lanthanum
      isotope_specification(58, 122, 175), & !cerium
      isotope_specification(59, 126, 178), & !praseodymium
      isotope_specification(60, 127, 178), & !neodymium
      isotope_specification(61, 130, 185), & !promethium
      isotope_specification(62, 133, 188), & !samarium
      isotope_specification(63, 136, 190), & !europium
      isotope_specification(64, 139, 191), & !gadolinium
      isotope_specification(65, 142, 192), & !terbium
      isotope_specification(66, 143, 193), & !dysprosium
      isotope_specification(67, 146, 196), & !holmium
      isotope_specification(68, 148, 196), & !erbium
      isotope_specification(69, 150, 198), & !thulium
      isotope_specification(70, 152, 200), & !ytterbium
      isotope_specification(71, 156, 209), & !lutetium
      isotope_specification(72, 159, 212), & !hafnium
      isotope_specification(73, 161, 217), & !tantalum
      isotope_specification(74, 163, 220), & !tungsten
      isotope_specification(75, 167, 225), & !rhenium
      isotope_specification(76, 169, 226), & !osmium
      isotope_specification(77, 172, 230), & !iridium
      isotope_specification(78, 175, 232), & !platinum
      isotope_specification(79, 178, 236), & !gold
      isotope_specification(80, 178, 239), & !mercury
      isotope_specification(81, 182, 245), & !thallium
      isotope_specification(82, 185, 246), & !lead
      isotope_specification(83, 188, 251), & !bismuth
      isotope_specification(84, 193, 237), & !polonium
      isotope_specification(85, 210, 211) & !astatine    
   ] 
    

    ! ------------------------------------------------------------------
    ! Phase 1: count total isotopes so we can allocate exact-sized arrays.
    ! We also need to know how many unique (Z, symbol) pairs there are so
    ! the element table can be sized.  Because Be appears in two groups
    ! (singleton Be7 + range Be9-10) we must deduplicate by Z.
    ! ------------------------------------------------------------------
    integer :: i, j, A, idx, eidx
    integer :: total_iso
    integer :: n_ele_unique
    logical :: found

    ! Count total isotopes across all spec groups
    total_iso = 0
    do i = 1, NISO
      total_iso = total_iso + (isotope_specs(i)%Amax - isotope_specs(i)%Amin + 1)
    end do

    ! Collect unique Z values for the element table (preserving first-seen order)
    n_ele_unique = 0
    do i = 1, NISO
      if (n_ele_unique < isotope_specs(i)%Z) then
       n_ele_unique = isotope_specs(i)%Z + 1 ! +1 for neutrons
      end if
    end do

    ! ------------------------------------------------------------------
    ! Phase 2: allocate and populate the tables
    ! ------------------------------------------------------------------
    if (allocated(isotope_table)) deallocate(isotope_table)
    if (allocated(element_table)) deallocate(element_table)

    allocate(isotope_table(total_iso))
    allocate(element_table(n_ele_unique))

    n_isotopes = total_iso
    n_elements = n_ele_unique

    ! Build element table first so we can look up indices cheaply below
    do j = 1, n_ele_unique
      element_table(j)%element = find_element_by_Z(j)
      element_table(j)%mass   = 0.0d0
    end do

    ! Build isotope table: iterate spec groups in order
    idx = 0
    do i = 1, NISO
      do A = isotope_specs(i)%Amin, isotope_specs(i)%Amax
        idx = idx + 1
        isotope_table(idx)%A      = A
        isotope_table(idx)%element = find_element_by_Z(isotope_specs(i)%Z)
        isotope_table(idx)%id     = make_iso_id(isotope_table(idx)%element, A)
        isotope_table(idx)%mass   = 0.0d0
      end do
    end do

    ! Cross-check: verify idx == total_iso
    if (idx /= total_iso) then
      write(0,'(a,i0,a,i0)') &
        'mesa_iso_mass WARNING: isotope count mismatch: built ', idx, &
        ', expected ', total_iso
    end if

    tables_built = .true.

  end subroutine init_iso_mass


  ! ===========================================================================
  ! reset_iso_mass
  !
  ! Zeroes all accumulated masses in both tables without rebuilding the
  ! isotope / element metadata.  Useful when processing multiple models
  ! sequentially without re-calling init_iso_mass.
  ! ===========================================================================
  subroutine reset_iso_mass()
    integer :: i

    if (.not. tables_built) then
      write(0,'(a)') 'mesa_iso_mass WARNING: reset called before init; calling init first.'
      call init_iso_mass()
      return
    end if

    do i = 1, n_isotopes
      isotope_table(i)%mass = 0.0d0
    end do
    do i = 1, n_elements
      element_table(i)%mass = 0.0d0
    end do
  end subroutine reset_iso_mass


  ! ===========================================================================
  ! add_iso_mass
  !
  ! Adds dm [g] to the accumulated mass of the isotope whose ID matches 'id'
  ! and propagates the same delta to the parent element.
  !
  ! Arguments:
  !   id (IN) : isotope ID string, e.g. "Fe56", "H1", "N1" (free neutron)
  !   dm (IN) : mass to add [g] (double precision)
  !
  ! Returns without error if 'id' is not found (a warning is written to
  ! stderr instead), so a missing species does not abort the calling loop.
  !
  ! Implementation:  binary search on isotope_table (sorted by id, which is
  ! lexicographically ordered because the symbol is title-cased and the mass
  ! number is zero-padded to the same width in the comparison key).
  ! Note: the table is ordered by network insertion order, not
  ! lexicographically, so we fall back to a linear search here.  For 3335
  ! isotopes this is still fast enough for typical zone counts (< 10⁵).
  ! ===========================================================================
  subroutine add_iso_mass(id, dm)
    character(len=*), intent(in) :: id
    real(real64),     intent(in) :: dm

    integer :: i, eidx
    logical :: found

    if (.not. tables_built) then
      write(0,'(a)') 'mesa_iso_mass ERROR: add_iso_mass called before init_iso_mass.'
      return
    end if

    found = .false.
    do i = 1, n_isotopes
      if (trim(isotope_table(i)%id) == trim(id)) then
        isotope_table(i)%mass = isotope_table(i)%mass + dm

        ! Propagate to parent element
        element_table(isotope_table(i)%element%Z + 1)%mass = element_table(isotope_table(i)%element%Z + 1)%mass + dm

        found = .true.
        exit
      end if
    end do

    if (.not. found) then
      write(0,'(3a)') 'mesa_iso_mass WARNING: isotope "', trim(id), &
                       '" not found in table; mass not added.'
    end if
  end subroutine add_iso_mass


  ! ===========================================================================
  ! get_iso_mass
  !
  ! Returns the currently accumulated mass [g] for isotope 'id'.
  ! Returns 0.0 and writes a warning if 'id' is not found.
  ! ===========================================================================
  function get_iso_mass(id) result(mass)
    character(len=*), intent(in) :: id
    real(real64) :: mass
    integer :: i

    mass = 0.0d0
    do i = 1, n_isotopes
      if (trim(isotope_table(i)%id) == trim(id)) then
        mass = isotope_table(i)%mass
        return
      end if
    end do
    write(0,'(3a)') 'mesa_iso_mass WARNING: get_iso_mass: "', trim(id), '" not found.'
  end function get_iso_mass

  ! ===========================================================================
  ! get_iso_mass
  !
  ! Returns the currently accumulated mass [g] for isotope 'id'.
  ! Returns 0.0 and writes a warning if 'id' is not found.
  ! ===========================================================================
  function get_iso_mass_by_idx(idx) result(mass)
    integer, intent(in) :: idx
    real(real64) :: mass
    integer :: i

    mass = 0.0d0
    if ((idx > 0).and.(idx <= n_isotopes)) then
        mass = isotope_table(idx)%mass
    else
      write(0,*) 'mesa_iso_mass WARNING: get_iso_mass_by_idx: "', idx, '" not found.'
    end if
    return
  end function get_iso_mass_by_idx


  ! ===========================================================================
  ! get_ele_mass
  !
  ! Returns the currently accumulated mass [g] for element with symbol 'sym'
  ! (title-cased, e.g. "Fe", "He").  Returns 0.0 and warns if not found.
  ! ===========================================================================
  function get_ele_mass(sym) result(mass)
    character(len=*), intent(in) :: sym
    real(real64) :: mass
    integer :: i

    mass = 0.0d0
    do i = 1, n_elements
      if (trim(element_table(i)%element%symbol) == trim(sym)) then
        mass = element_table(i)%mass
        return
      end if
    end do
    write(0,'(3a)') 'mesa_iso_mass WARNING: get_ele_mass: "', trim(sym), '" not found.'
  end function get_ele_mass


  ! ===========================================================================
  ! get_ele_mass
  !
  ! Returns the currently accumulated mass [g] for element with symbol 'sym'
  ! (title-cased, e.g. "Fe", "He").  Returns 0.0 and warns if not found.
  ! ===========================================================================
  function get_ele_mass_by_Z(Z) result(mass)
    integer, intent(in) :: Z
    real(real64) :: mass

    mass = 0.0d0
    if (Z < n_elements) then
        mass = element_table(Z + 1)%mass
    else
      write(0,'(3a)') 'mesa_iso_mass WARNING: get_ele_mass_by_Z: Z = ', Z, '" not found.'
    end if
  end function get_ele_mass_by_Z

  ! ===========================================================================
  ! PRIVATE HELPERS
  ! ===========================================================================

  ! ---------------------------------------------------------------------------
  ! make_iso_id
  !
  ! Constructs the title-cased isotope ID string used by mesa_mod_summary,
  ! e.g. ("fe", 56) -> "Fe56", ("n", 1) -> "N1".
  !
  ! The free neutron ("n", Z=0, A=1) is a special case: MESA's column header
  ! is "neut" but mesa_mod_summary already normalises it to "N1", which is
  ! what we produce here.
  ! ---------------------------------------------------------------------------
  function make_iso_id(elem, A) result(id)
    type(element_type), intent(in) :: elem
    integer,          intent(in) :: A
    character(len=6)             :: id

    character(len=3) :: tc
    character(len=4) :: anum

    write(anum, '(i0)') A
    id = trim(elem%symbol) // trim(anum)
  end function make_iso_id


  ! ---------------------------------------------------------------------------
  ! title_case
  !
  ! Returns a copy of 'sym' with the first character uppercased and the rest
  ! unchanged.  Input is expected to be a lower-case MESA symbol (1-3 chars).
  ! ---------------------------------------------------------------------------
  function title_case(sym) result(tc)
    character(len=*), intent(in) :: sym
    character(len=len(sym))      :: tc
    integer :: offset

    tc = adjustl(sym)

    ! Convert first character to upper case (ASCII offset 'a'-'A' = 32)
    offset = iachar(tc(1:1)) - iachar('a')
    if (offset >= 0 .and. offset <= 25) then
      tc(1:1) = achar(iachar(tc(1:1)) - 32)
    end if
  end function title_case


  ! ---------------------------------------------------------------------------
  ! first_spec_for_Z
  !
  ! Returns the index of the first entry in ele_Z(:) whose value equals
  ! target_Z.  Used during init to find the symbol for each unique element.
  ! ---------------------------------------------------------------------------
  function first_spec_for_Z(target_Z, ele_Z, n) result(idx)
    integer, intent(in) :: target_Z
    integer, intent(in) :: n
    integer, intent(in) :: ele_Z(n)
    integer             :: idx, i

    idx = 1
    do i = 1, n
      if (ele_Z(i) == target_Z) then
        idx = i
        return
      end if
    end do
  end function first_spec_for_Z

end module mesa_iso_mass
