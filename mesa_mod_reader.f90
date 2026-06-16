! ==============================================================================
! MODULE: mesa_mod_reader
!
! Parses a MESA stellar model file (.mod) and exposes its contents through the
! mesa_model derived type.
!
! MESA .mod file structure
! ------------------------
!   Line  1-2  : Comment lines beginning with '!'
!   Line  3    : 12-space prefix, integer model version, then descriptive text
!   Line  4    : Blank
!   Lines 5-?  : Pre-table key-value metadata (label in chars 1-33, value
!                after that); ends at the first blank line.  Typical entries
!                include "M/Msun", "R/Rsun", "L/Lsun", "species", "n_shells".
!   Next line  : Column header row (whitespace-delimited token strings)
!   Data rows  : One row per stellar zone, outermost first.
!                  col 1        : 5-char zone number
!                  cols 2..ncols: fixed-width 27-char fields in D-exponent
!                                 engineering notation
!                Terminated by a blank line.
!   Next line  : Literal "        previous model" sentinel
!   Blank line
!   4 kv lines : Post-table metadata in the same label/value format
!   Blank / EOF
!
! Public interface
! ----------------
!   mesa_model         – derived type holding the entire parsed model
!   read_mesa_model    – subroutine: filename → mesa_model
!   get_meta           – subroutine: look up a metadata key by label string
!   destroy_mesa_model – subroutine: frees all allocatable components
!
! Notes
! -----
!   * kv_pair values are stored as raw strings so that numeric conversions
!     are left to the caller; this insulates the reader from future format
!     changes.
!   * The data table is stored as real(kind=8) (double precision).
!   * Row 1 of model%table is the outermost stellar zone; row nrows is the
!     stellar centre.
! ==============================================================================
module mesa_mod_reader
  use iso_fortran_env, only: real64
  implicit none
  private

  ! ---------------------------------------------------------------------------
  ! Module-level format constants
  ! These mirror the fixed-width column layout used in MESA .mod files.
  ! ---------------------------------------------------------------------------
  integer, parameter, public :: MESA_LABEL_LEN   = 33   ! fixed width of metadata label field
  integer, parameter, public :: MESA_VALUE_LEN   = 200  ! max length of a metadata value string
  integer, parameter, public :: MESA_HEADER_LEN  = 33   ! max length of a column header token
  integer, parameter, public :: MESA_COL1_WIDTH  = 5    ! width of the leading zone-number column
  integer, parameter, public :: MESA_DATACOL_WIDTH = 27 ! width of each numeric data column

  ! ---------------------------------------------------------------------------
  ! kv_pair: one key-value metadata entry
  !
  ! Both fields are stored as trimmed strings.  The reader strips surrounding
  ! single quotes from quoted string values so callers can do a direct
  ! internal read without having to strip them again.
  ! ---------------------------------------------------------------------------
  type, public :: kv_pair
    character(len=MESA_LABEL_LEN) :: label = ''   ! metadata key  (e.g. "M/Msun")
    character(len=MESA_VALUE_LEN) :: value = ''   ! metadata value (e.g. "0.8500000")
  end type kv_pair

  ! ---------------------------------------------------------------------------
  ! mesa_model: complete in-memory representation of one MESA .mod file
  !
  ! After a successful call to read_mesa_model, all allocatable components
  ! are populated.  ncols == -1 signals a parse error.
  !
  ! Table layout  (model%table(irow, icol)):
  !   icol 1  : ln(rho)  – natural log of density       [ln g/cm³]
  !   icol 2  : ln(T)    – natural log of temperature   [ln K]
  !   icol 3  : ln(R)    – natural log of radius        [ln cm]
  !   icol 4  : L        – luminosity
  !   icol 5  : v        – velocity
  !   icol 6  : (additional structural quantity)
  !   icol 7+ : mass fractions for each isotope in the nuclear network
  !             (order matches model%col_header(7:ncols))
  ! ---------------------------------------------------------------------------
  type, public :: mesa_model

    integer :: mesa_model_version = 0   ! format version integer from line 3

    ! Pre-table metadata (label/value pairs from lines 5 onward)
    integer                    :: n_meta          = 0   ! number of kv pairs stored
    integer                    :: meta_alloc_size = 0   ! current allocated capacity
    type(kv_pair), allocatable :: meta(:)               ! size n_meta after finalisation

    ! Data table dimensions
    integer :: nrows        = 0   ! number of stellar zones
    integer :: ncols        = 0   ! number of data columns  (-1 on error)
    integer :: num_isotopes = 0   ! number of isotope columns  (= ncols - 6)

    integer,          allocatable :: row_number(:)    ! MESA zone numbers (size nrows)
    character(len=MESA_HEADER_LEN), &
                      allocatable :: col_header(:)    ! column header strings (size ncols)
    real(kind=real64),allocatable :: table(:,:)       ! data values (nrows × ncols)

    ! Post-table metadata (4 kv pairs after the "previous model" sentinel)
    integer                    :: n_post_meta          = 0
    integer                    :: post_meta_alloc_size = 0
    type(kv_pair), allocatable :: post_meta(:)

  end type mesa_model

  ! ---------------------------------------------------------------------------
  ! Public procedure interfaces
  ! ---------------------------------------------------------------------------
  public :: read_mesa_model
  public :: get_meta
  public :: destroy_mesa_model


contains

  ! ===========================================================================
  ! SECTION: Dynamic metadata array management
  !
  ! The pre-table metadata block has a variable number of entries (typically
  ! ~18), and the post-table block has exactly 4.  Both use a simple doubling
  ! strategy so that parse_kv_line can be called in a read-until-blank loop
  ! without knowing the count in advance.
  !
  ! Each group has three helpers:
  !   add_kv / add_kv_post       – append one kv_pair, reallocating if full
  !   realloc_meta_kv / ...      – double the backing array capacity
  !   finalize_meta_size / ...   – shrink the backing array to exact count
  ! ===========================================================================

  ! ---------------------------------------------------------------------------
  ! add_kv: append one kv_pair to model%meta, growing the array if needed
  ! ---------------------------------------------------------------------------
  subroutine add_kv(model, pair)
    type(mesa_model), intent(inout) :: model
    type(kv_pair),    intent(in)    :: pair

    if (model%meta_alloc_size == 0 .or. &
        model%meta_alloc_size == model%n_meta) then
      call realloc_meta_kv(model)
    end if
    model%n_meta = model%n_meta + 1
    model%meta(model%n_meta) = pair
  end subroutine add_kv

  ! ---------------------------------------------------------------------------
  ! realloc_meta_kv: double the capacity of model%meta
  !   First call: allocates 20 slots.
  !   Subsequent calls: copies existing data into a 2× larger array.
  ! ---------------------------------------------------------------------------
  subroutine realloc_meta_kv(model)
    type(mesa_model), intent(inout) :: model
    type(kv_pair), allocatable :: kv_pairs(:)
    integer :: ii, newAllocSize

    if (model%meta_alloc_size == 0) then
      allocate(model%meta(20))
      model%meta_alloc_size = 20
    else
      newAllocSize = model%meta_alloc_size * 2
      allocate(kv_pairs(newAllocSize))
      do ii = 1, model%n_meta
        kv_pairs(ii) = model%meta(ii)
      end do
      deallocate(model%meta)
      model%meta = kv_pairs
      model%meta_alloc_size = newAllocSize
    end if
  end subroutine realloc_meta_kv

  ! ---------------------------------------------------------------------------
  ! finalize_meta_size: shrink model%meta to exactly n_meta elements
  !   Called once after all pre-table metadata has been read so that the
  !   array does not carry unused trailing slots.
  ! ---------------------------------------------------------------------------
  subroutine finalize_meta_size(model)
    type(mesa_model), intent(inout) :: model
    type(kv_pair), allocatable :: kv_pairs(:)
    integer :: ii

    allocate(kv_pairs(model%n_meta))
    do ii = 1, model%n_meta
      kv_pairs(ii) = model%meta(ii)
    end do
    deallocate(model%meta)
    model%meta = kv_pairs
    model%meta_alloc_size = model%n_meta
  end subroutine finalize_meta_size

  ! ---------------------------------------------------------------------------
  ! add_kv_post: append one kv_pair to model%post_meta (same logic as add_kv)
  ! ---------------------------------------------------------------------------
  subroutine add_kv_post(model, pair)
    type(mesa_model), intent(inout) :: model
    type(kv_pair),    intent(in)    :: pair

    if (model%post_meta_alloc_size == 0 .or. &
        model%post_meta_alloc_size == model%n_post_meta) then
      call realloc_post_meta_kv(model)
    end if
    model%n_post_meta = model%n_post_meta + 1
    model%post_meta(model%n_post_meta) = pair
  end subroutine add_kv_post

  ! ---------------------------------------------------------------------------
  ! realloc_post_meta_kv: double capacity of model%post_meta
  !   First call: allocates 4 slots (enough for the typical 4-line block).
  ! ---------------------------------------------------------------------------
  subroutine realloc_post_meta_kv(model)
    type(mesa_model), intent(inout) :: model
    type(kv_pair), allocatable :: kv_pairs(:)
    integer :: ii, newAllocSize

    if (model%post_meta_alloc_size == 0) then
      allocate(model%post_meta(4))
      model%post_meta_alloc_size = 4
    else
      newAllocSize = model%post_meta_alloc_size * 2
      allocate(kv_pairs(newAllocSize))
      do ii = 1, model%n_post_meta
        kv_pairs(ii) = model%post_meta(ii)
      end do
      deallocate(model%post_meta)
      model%post_meta = kv_pairs
      model%post_meta_alloc_size = newAllocSize
    end if
  end subroutine realloc_post_meta_kv

  ! ---------------------------------------------------------------------------
  ! finalize_post_meta_size: shrink model%post_meta to exactly n_post_meta
  ! ---------------------------------------------------------------------------
  subroutine finalize_post_meta_size(model)
    type(mesa_model), intent(inout) :: model
    type(kv_pair), allocatable :: kv_pairs(:)
    integer :: ii

    allocate(kv_pairs(model%n_post_meta))
    do ii = 1, model%n_post_meta
      kv_pairs(ii) = model%post_meta(ii)
    end do
    deallocate(model%post_meta)
    model%post_meta = kv_pairs
    model%post_meta_alloc_size = model%n_post_meta
  end subroutine finalize_post_meta_size


  ! ===========================================================================
  ! get_meta
  !
  ! Searches model%meta for the first entry whose label matches 'id' and
  ! returns it in kv_pair_result.  If no match is found, kv_pair_result
  ! is left with its default-initialised (empty) values.
  !
  ! Arguments:
  !   model          (IN)  : parsed MESA model
  !   id             (IN)  : label to search for (e.g. "M/Msun")
  !   kv_pair_result (OUT) : matching kv_pair, or empty if not found
  ! ===========================================================================
  subroutine get_meta(model, id, kv_pair_result)
    type(mesa_model), intent(in)  :: model
    character(len=*), intent(in)  :: id
    type(kv_pair),    intent(out) :: kv_pair_result

    integer :: i

    do i = 1, model%n_meta
      if (model%meta(i)%label == id) then
        kv_pair_result = model%meta(i)
        return   ! stop at first match
      end if
    end do
  end subroutine get_meta


  ! ===========================================================================
  ! read_mesa_model
  !
  ! Opens 'filename', parses it completely following the MESA .mod layout
  ! described in the module header, and returns the result in 'model'.
  !
  ! On any I/O or parse error a message is written to stderr and
  ! model%ncols is set to -1 before returning.
  !
  ! Arguments:
  !   filename (IN)  : path to the .mod file
  !   model    (OUT) : populated mesa_model on success; ncols == -1 on error
  ! ===========================================================================
  subroutine read_mesa_model(filename, model)
    character(len=*), intent(in)  :: filename
    type(mesa_model), intent(out) :: model

    integer, parameter :: UNIT     = 42      ! logical unit number for file I/O
    integer, parameter :: MAX_ROWS = 200000  ! safety cap on zone count
    integer, parameter :: N_POST   = 4       ! number of post-table metadata lines

    ! Allocatable line buffer: resized after species count is known so it can
    ! hold one full data row without truncation.
    character(len=:), allocatable :: line

    integer        :: ios, i, icol, irow
    integer        :: mesa_model_version
    integer        :: ncols
    integer        :: num_species   ! isotope count read from "species" metadata key
    integer        :: num_data      ! zone count read from "n_shells" metadata key
    type(kv_pair)  :: pair1, pair2
    logical        :: end_of_meta
    logical        :: haspair2
    character(len=50) :: tmpstring

    ! Temporary storage for table rows; exact-sized arrays are allocated after
    ! nrows and ncols are known, then populated in a second pass.
    integer,       allocatable :: tmp_rownum(:)
    real(real64),  allocatable :: tmp_table(:,:)

    ! ------------------------------------------------------------------
    ! Open file
    ! ------------------------------------------------------------------
    ! Initial line buffer large enough for the header lines (which can be wide
    ! due to metadata values) before we know the data column count.
    allocate(character(len=50000) :: line)

    open(unit=UNIT, file=trim(filename), status='old', action='read', &
         iostat=ios)
    if (ios /= 0) then
      write(0,'(3a)') 'mesa_mod_reader ERROR: cannot open "', &
                       trim(filename), '"'
      model%ncols = -1
      return
    end if

    ! ------------------------------------------------------------------
    ! Lines 1-2: comment lines (expected to begin with '!')
    ! ------------------------------------------------------------------
    do i = 1, 2
      read(UNIT, '(a)', iostat=ios) line
      if (ios /= 0) then
        call io_error('reading comment line', i)
        model%ncols = -1; return
      end if
      if (len_trim(line) > 0 .and. line(1:1) /= '!') &
        write(0,'(a,i0)') 'mesa_mod_reader WARNING: expected ! on line ', i
    end do

    ! ------------------------------------------------------------------
    ! Line 3: model version integer
    !   Format: 12 leading spaces, integer, then descriptive text.
    !   We skip the first 12 characters and read the integer from what
    !   remains so that the text suffix does not cause a parse error.
    ! ------------------------------------------------------------------
    read(UNIT, '(a)', iostat=ios) line
    if (ios /= 0) then
      call io_error('reading model version line', 3)
      model%ncols = -1; return
    end if
    read(line(13:), *, iostat=ios) mesa_model_version
    if (ios /= 0) then
      write(0,'(a)') 'mesa_mod_reader ERROR: cannot parse model version from line 3'
      model%mesa_model_version = -1; return
    end if
    model%mesa_model_version = mesa_model_version

    ! ------------------------------------------------------------------
    ! Line 4: blank separator — read and discard
    ! ------------------------------------------------------------------
    read(UNIT, '(a)', iostat=ios) line

    ! ------------------------------------------------------------------
    ! Lines 5+: pre-table key-value metadata
    !
    ! We read in a loop until we hit a blank line (end_of_meta = .true.).
    ! Each non-blank line is parsed by parse_kv_line into up to two
    ! kv_pair values and appended to the dynamic model%meta array.
    !
    ! Two keys have special significance to the reader itself:
    !   "species"  -> num_species : number of isotopes in the network
    !   "n_shells" -> num_data    : number of stellar zones (data rows)
    ! Both are needed to size the line buffer and allocate the table.
    ! ------------------------------------------------------------------
    end_of_meta = .false.
    num_species = 0
    num_data    = 0

    do while (.not. end_of_meta)
      read(UNIT, '(a)', iostat=ios) line
      if (ios /= 0) then
        call io_error('reading metadata line', -1)
        model%ncols = -1; return
      end if

      end_of_meta = (len(trim(line)) == 0)

      if (.not. end_of_meta) then
        call parse_kv_line(line, pair1, haspair2, pair2)
        call add_kv(model, pair1)

        ! Capture species count for later table allocation
        if (pair1%label == "species") then
          tmpstring = trim(pair1%value)
          read(tmpstring, *) num_species
        end if

        ! Capture zone count for later table allocation
        if (pair1%label == "n_shells") then
          tmpstring = trim(pair1%value)
          read(tmpstring, *) num_data
        end if

        ! Some metadata lines carry a second label/value pair inline
        if (haspair2) call add_kv(model, pair2)
      end if
    end do

    call finalize_meta_size(model)
    model%num_isotopes = num_species
    model%nrows        = num_data

    ! Resize line buffer: each data row is  MESA_COL1_WIDTH + ncols*MESA_DATACOL_WIDTH
    ! plus a small safety margin.
    deallocate(line)
    allocate(character(len=(num_species + 6) * MESA_DATACOL_WIDTH &
                            + MESA_COL1_WIDTH + 10) :: line)

    ! ------------------------------------------------------------------
    ! Column header line
    !   Immediately follows the blank line that ended the metadata block.
    !   Contains whitespace-delimited token strings, one per data column.
    ! ------------------------------------------------------------------
    read(UNIT, '(a)', iostat=ios) line
    if (ios /= 0) then
      call io_error('reading column header line', -1)
      model%ncols = -1; return
    end if

    ! ncols = 6 fixed structural columns + one column per isotope species
    ncols        = num_species + 6
    model%ncols  = ncols
    allocate(model%col_header(ncols))
    call parse_column_headers(line, ncols, model%col_header)

    ! ------------------------------------------------------------------
    ! Data table
    !
    ! Each row represents one stellar zone; zones are listed outermost
    ! first.  The row begins with a 5-character zone number followed by
    ! ncols fixed-width (27-char) fields in D-exponent notation.
    !
    ! We preallocate temporary arrays sized to num_data (from metadata)
    ! and fill them row by row, stopping at the first blank line.
    ! The exact row count is then stored in model%nrows and exact-sized
    ! final arrays are allocated and filled from the temporaries.
    ! ------------------------------------------------------------------
    allocate(tmp_rownum(num_data))
    allocate(tmp_table(num_data, ncols))

    irow = 0
    do
      read(UNIT, '(a)', iostat=ios) line
      if (ios /= 0)            exit   ! EOF or read error
      if (len_trim(line) == 0) exit   ! blank line marks end of data table

      irow = irow + 1
      if (irow > MAX_ROWS) then
        write(0,'(a)') 'mesa_mod_reader ERROR: row count exceeds MAX_ROWS'
        model%ncols = -1; return
      end if

      ! Parse the leading zone-number field
      read(line(1:MESA_COL1_WIDTH), *, iostat=ios) tmp_rownum(irow)
      if (ios /= 0) then
        write(0,'(a,i0)') &
          'mesa_mod_reader WARNING: bad row number at table row ', irow
        tmp_rownum(irow) = irow   ! fall back to sequential numbering
      end if

      ! Parse each numeric data column
      do icol = 1, ncols
        call read_table_value(line, MESA_COL1_WIDTH, MESA_DATACOL_WIDTH, &
                              icol, tmp_table(irow, icol), ios)
        if (ios /= 0) then
          write(0,'(a,i0,a,i0)') &
            'mesa_mod_reader WARNING: parse error at row ', irow, ' col ', icol
          tmp_table(irow, icol) = 0.0d0
        end if
      end do
    end do

    ! Copy exact-sized slices into the model (overwrite nrows with actual count)
    model%nrows = irow
    allocate(model%row_number(irow))
    allocate(model%table(irow, ncols))
    model%row_number(1:irow)         = tmp_rownum(1:irow)
    model%table(1:irow, 1:ncols)     = tmp_table(1:irow, 1:ncols)
    deallocate(tmp_rownum, tmp_table)

    ! ------------------------------------------------------------------
    ! Post-table section
    !
    ! Layout after the data-table blank line:
    !   "        previous model"   (sentinel line, consumed above by the
    !                               blank-line exit from the data loop)
    !   Blank line
    !   4 (or more) key-value metadata lines in the same format as the
    !   pre-table block.
    !   Blank line / EOF
    !
    ! The "previous model" sentinel was already consumed as the blank-
    ! line terminator of the data loop, so we read two more lines here
    ! (the sentinel itself and its trailing blank) before the kv block.
    ! ------------------------------------------------------------------

    ! Read the "        previous model" sentinel line
    read(UNIT, '(a)', iostat=ios) line
    if (ios /= 0) then
      call io_error('reading "previous model" sentinel', -1)
      model%ncols = -1; return
    end if

    ! Blank line after the sentinel
    read(UNIT, '(a)', iostat=ios) line

    ! Read post-table key-value pairs until another blank line or EOF
    end_of_meta = .false.
    do while (.not. end_of_meta)
      read(UNIT, '(a)', iostat=ios) line
      if (ios /= 0) exit   ! EOF is a valid terminator here

      end_of_meta = (len(trim(line)) == 0)
      if (.not. end_of_meta) then
        call parse_kv_line(line, pair1, haspair2, pair2)
        call add_kv_post(model, pair1)
      end if
    end do
    call finalize_post_meta_size(model)

    close(UNIT)

  end subroutine read_mesa_model


  ! ===========================================================================
  ! destroy_mesa_model
  !
  ! Releases all allocatable components of a mesa_model and resets all
  ! integer counters to zero.  Call this when the model is no longer needed
  ! to avoid memory leaks.
  ! ===========================================================================
  subroutine destroy_mesa_model(model)
    type(mesa_model), intent(inout) :: model

    if (allocated(model%meta))       deallocate(model%meta)
    if (allocated(model%col_header)) deallocate(model%col_header)
    if (allocated(model%row_number)) deallocate(model%row_number)
    if (allocated(model%table))      deallocate(model%table)
    if (allocated(model%post_meta))  deallocate(model%post_meta)

    model%ncols                = 0
    model%n_meta               = 0
    model%nrows                = 0
    model%n_post_meta          = 0
    model%meta_alloc_size      = 0
    model%post_meta_alloc_size = 0

  end subroutine destroy_mesa_model


  ! ===========================================================================
  ! PRIVATE HELPERS
  ! ===========================================================================

  ! ---------------------------------------------------------------------------
  ! parse_kv_line
  !
  ! Parses one metadata line into up to two key-value pairs.
  !
  ! Line format:
  !   <MESA_LABEL_LEN chars : label> <whitespace> <value token>
  !   [<label2 token> <whitespace> <value2 token>]
  !
  ! The label occupies the first MESA_LABEL_LEN characters (leading/trailing
  ! spaces stripped).  The value follows and is read by extract_value, which
  ! handles quoted strings and D/E-exponent numeric tokens.
  !
  ! Arguments:
  !   line     (IN)  : raw line from the file
  !   pair     (OUT) : primary key-value pair
  !   hasPair2 (OUT) : .true. if a second pair was found
  !   pair2    (OUT) : secondary key-value pair (valid only when hasPair2)
  ! ---------------------------------------------------------------------------
  subroutine parse_kv_line(line, pair, hasPair2, pair2)
    character(len=*), intent(in)  :: line
    type(kv_pair),    intent(out) :: pair
    logical,          intent(out) :: hasPair2
    type(kv_pair),    intent(out) :: pair2

    integer :: pos, llen, ii, endidx
    character(len=len(line)) :: rest

    ! Initialise outputs
    pair%label  = '';  pair%value  = ''
    pair2%label = '';  pair2%value = ''
    hasPair2 = .false.

    llen = len_trim(line)
    if (llen < 1) return

    ! Extract primary label from the first MESA_LABEL_LEN characters
    pair%label = adjustl(trim(line(1 : min(MESA_LABEL_LEN, llen))))
    if (llen <= MESA_LABEL_LEN) return   ! line ends before any value

    ! Extract primary value from the remainder of the line
    rest = adjustl(line(MESA_LABEL_LEN + 1 :))
    call extract_value(rest, pair%value, pos)

    ! Check whether a second label+value pair follows on the same line
    if (pos <= len_trim(rest)) then
      rest = adjustl(trim(rest(pos:)))
      if (len(rest) > 1) then
        ! Find the end of the second label (next space)
        endidx = len_trim(rest)
        do ii = 1, len(rest)
          if (iachar(rest(ii:ii)) == 32) then   ! space character
            endidx = ii - 1
            exit
          end if
        end do
        hasPair2    = .true.
        pair2%label = rest(1 : endidx)
        rest        = adjustl(rest(endidx + 1:))
        call extract_value(rest, pair2%value, pos)
      end if
    end if

  end subroutine parse_kv_line


  ! ---------------------------------------------------------------------------
  ! extract_value
  !
  ! Extracts the first value token from 'str' (after stripping leading spaces)
  ! and returns the position of the first character after the token in
  ! 'next_pos' (1-based, relative to the leading-space-stripped copy).
  !
  ! Token types handled:
  !   'quoted string' – everything between the first pair of single quotes,
  !                     quotes themselves are not included in val
  !   unquoted token  – characters up to the next space or end of string
  !
  ! Arguments:
  !   str      (IN)  : input string (may have leading spaces)
  !   val      (OUT) : extracted value token (quotes stripped if applicable)
  !   next_pos (OUT) : index of the next character after the token in str
  ! ---------------------------------------------------------------------------
  subroutine extract_value(str, val, next_pos)
    character(len=*), intent(in)  :: str
    character(len=*), intent(out) :: val
    integer,          intent(out) :: next_pos

    integer :: qend, spc
    character(len=len(str)) :: s

    val      = ''
    next_pos = len_trim(str) + 1

    s = adjustl(str)
    if (len_trim(s) == 0) return

    if (s(1:1) == "'") then
      ! Quoted string: find the matching closing quote
      qend = index(s(2:), "'")
      if (qend == 0) then
        ! No closing quote — consume the rest of the string
        val      = trim(s(2:))
        next_pos = len_trim(s) + 1
      else
        val      = s(2 : qend)   ! content between the quotes (qend is 1-based in s(2:))
        next_pos = qend + 2      ! position after the closing quote in s
      end if
    else
      ! Unquoted token: ends at the next space or end of string
      spc = index(trim(s), ' ')
      if (spc == 0) then
        val      = trim(s)
        next_pos = len_trim(s) + 1
      else
        val      = s(1 : spc - 1)
        next_pos = spc + 1
      end if
    end if

  end subroutine extract_value


  ! ---------------------------------------------------------------------------
  ! parse_column_headers
  !
  ! Splits 'line' into up to n whitespace-delimited tokens and stores them
  ! in headers(1:n).  Tokens shorter than MESA_HEADER_LEN are right-padded
  ! with spaces by the assignment to the fixed-length character array.
  !
  ! Arguments:
  !   line    (IN)  : the column-header line from the .mod file
  !   n       (IN)  : expected number of column headers (= model%ncols)
  !   headers (OUT) : array of parsed header strings
  ! ---------------------------------------------------------------------------
  subroutine parse_column_headers(line, n, headers)
    character(len=*),               intent(in)  :: line
    integer,                        intent(in)  :: n
    character(len=MESA_HEADER_LEN), intent(out) :: headers(n)

    character(len=len(line)) :: rest
    integer :: icol, spc

    headers = ''
    rest    = adjustl(line)

    do icol = 1, n
      rest = adjustl(rest)
      if (len_trim(rest) == 0) exit   ! fewer tokens than expected

      spc = index(trim(rest), ' ')
      if (spc == 0) then
        headers(icol) = trim(rest)
        rest          = ''
      else
        headers(icol) = rest(1 : spc - 1)
        rest          = rest(spc + 1 :)
      end if
    end do

  end subroutine parse_column_headers


  ! ---------------------------------------------------------------------------
  ! read_table_value
  !
  ! Slices the icol-th data field from a table row line and parses it as a
  ! double-precision floating-point number.
  !
  ! Column layout within a data row:
  !   Characters 1 .. col1_w           : zone number (integer)
  !   Characters col1_w+1 .. col1_w+col_w  : column 1
  !   Characters col1_w+col_w+1 .. ...     : column 2  etc.
  !
  ! The MESA file uses D-exponent notation (e.g. 1.23D+00).  Fortran's
  ! list-directed internal read accepts D as well as E, so no substitution
  ! is needed.
  !
  ! Arguments:
  !   line   (IN)  : one raw data row from the file
  !   col1_w (IN)  : width of the leading zone-number column
  !   col_w  (IN)  : fixed width of each numeric data column
  !   icol   (IN)  : 1-based column index to extract
  !   val    (OUT) : parsed value
  !   ios    (OUT) : I/O status; non-zero indicates a parse failure
  ! ---------------------------------------------------------------------------
  subroutine read_table_value(line, col1_w, col_w, icol, val, ios)
    character(len=*), intent(in)  :: line
    integer,          intent(in)  :: col1_w, col_w, icol
    real(real64),     intent(out) :: val
    integer,          intent(out) :: ios

    integer :: cstart, cend
    character(len=col_w) :: token

    cstart = col1_w + (icol - 1) * col_w + 1
    cend   = cstart + col_w - 1

    if (cend > len(line)) then
      ! Line is shorter than expected — column is missing
      ios = -1
      val = 0.0d0
      return
    end if

    token = line(cstart : cend)
    read(token, *, iostat=ios) val

  end subroutine read_table_value


  ! ---------------------------------------------------------------------------
  ! io_error
  !
  ! Writes a formatted error message to stderr (unit 0).
  ! Pass lineno <= 0 to omit the line-number suffix.
  ! ---------------------------------------------------------------------------
  subroutine io_error(msg, lineno)
    character(len=*), intent(in) :: msg
    integer,          intent(in) :: lineno

    if (lineno > 0) then
      write(0,'(a,a,a,i0)') 'mesa_mod_reader ERROR: ', trim(msg), &
                              ' at file line ', lineno
    else
      write(0,'(a,a)') 'mesa_mod_reader ERROR: ', trim(msg)
    end if
  end subroutine io_error

end module mesa_mod_reader
