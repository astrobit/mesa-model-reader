! ==============================================================================
! MODULE: unique_strings_mod
!
! Provides a single utility subroutine for deduplicating a character array
! while preserving first-occurrence order.  Used by mesa_mod_summary to
! collapse the full isotope list down to its set of unique element symbols.
! ==============================================================================
module unique_strings_mod
  implicit none

contains

  !-----------------------------------------------------------------------------
  ! SUBROUTINE: find_unique_strings
  !
  ! Scans an array of fixed-length character strings and returns a pointer to
  ! a freshly allocated array containing only the first occurrence of each
  ! distinct value (order preserved).
  !
  ! Arguments:
  !   input      (IN)  : character array to scan
  !   n_input    (IN)  : number of elements in input
  !   output     (OUT) : pointer to allocated array of unique strings
  !   n_unique   (OUT) : number of unique strings found (= size of output)
  !
  ! Caller is responsible for deallocating output when finished.
  !-----------------------------------------------------------------------------
  subroutine find_unique_strings(input, n_input, output, n_unique)
    integer,                         intent(in)  :: n_input
    character(len=*),                intent(in)  :: input(n_input)
    character(len=len(input)), pointer, &
                                     intent(out) :: output(:)
    integer,                         intent(out) :: n_unique

    ! Temporary workspace sized for the worst case (all strings unique)
    character(len=len(input)), allocatable :: tmp(:)
    integer :: i, j
    logical :: is_duplicate

    allocate(tmp(n_input))
    n_unique = 0

    do i = 1, n_input
      is_duplicate = .false.

      ! Compare against every string already accepted as unique
      do j = 1, n_unique
        if (trim(input(i)) == trim(tmp(j))) then
          is_duplicate = .true.
          exit
        end if
      end do

      ! First occurrence: add to temporary buffer
      if (.not. is_duplicate) then
        n_unique = n_unique + 1
        tmp(n_unique) = input(i)
      end if
    end do

    ! Copy exactly the right number of results into the output pointer
    allocate(output(n_unique))
    output(1:n_unique) = tmp(1:n_unique)

  end subroutine find_unique_strings

end module unique_strings_mod

module helper_functions
  use iso_fortran_env, only: real64
  implicit none

contains

  subroutine printMassSolar(mass, strOut)
  real(kind=real64), intent(in) :: mass
  character(len=*), intent(out) :: strOut
  real, parameter :: solar_mass  = (1.32712442099e20 / 6.67430e-11) * 1.0e3
    if (log10(mass) < 30) then
      write(strOut, '(E14.7)') mass / solar_mass
    else
      write(strOut, '(F11.7)') mass / solar_mass
    end if
  end subroutine printMassSolar
end module helper_functions

! ==============================================================================
! PROGRAM: mesa_mod_summary
!
! Reads a MESA stellar model file (.mod) via mesa_mod_reader and prints a
! concise physical summary:
!
!   - Maximum and central log temperatures  [log K]
!   - Outer radius in cm, Earth radii, or solar radii
!   - Numerically integrated total mass     [g and M☉]
!   - Mass of each element present          [M☉ and % of total]
!   - Tab-separated data rows ready for comparison spreadsheets (CO WD and
!     nuclear-network output tables)
!
! Usage:
!   ./mesa_mod_summary  <path-to-model.mod>
!
! Table column layout (from mesa_mod_reader):
!   col 1 : ln(rho)  – log density      [ln g/cm³]
!   col 2 : ln(T)    – log temperature  [ln K]
!   col 3 : ln(R)    – log radius       [ln cm]
!   col 4 : L        – luminosity
!   col 5 : v        – velocity
!   col 6 : (unused in this program)
!   col 7+ : mass fractions for each isotope in the nuclear network
!
! Isotope naming convention used by MESA:
!   Single-character element symbol followed by mass number, e.g. "h1", "he4"
!   Special cases handled here: "neut" -> "n1", "prot" -> "H1"
! ==============================================================================
program mesa_mod_summary
  use unique_strings_mod
  use mesa_mod_reader
  use mesa_iso_mass
  use helper_functions
  use iso_fortran_env, only: real64

  implicit none

  ! ---------------------------------------------------------------------------
  ! Derived-type instance holding the entire parsed model
  ! ---------------------------------------------------------------------------
  type(mesa_model) :: model

  ! ---------------------------------------------------------------------------
  ! Command-line argument handling
  ! ---------------------------------------------------------------------------
  integer :: num_args, ix, jx, irow, icol, nm
  character(len=256), dimension(:), allocatable :: args

  ! ---------------------------------------------------------------------------
  ! Physical constants and unit conversion factors
  ! ---------------------------------------------------------------------------
  real, parameter :: PI          = 4.0 * atan(1.0)    ! π
  real, parameter :: solar_mass  = (1.32712442099e20 / 6.67430e-11) * 1.0e3
                                                       ! M☉ in grams
  real, parameter :: solar_radius = 6.957e10           ! R☉ in cm
  real, parameter :: earth_radius = 6.371e8            ! R⊕ in cm

  ! ---------------------------------------------------------------------------
  ! Scalar physical quantities derived from the model
  ! ---------------------------------------------------------------------------
  real(kind=real64)    :: radius           ! outer radius of the model [cm]
  real(kind=real64)    :: radius_rel_earth ! outer radius in Earth radii
  real(kind=real64)    :: radius_rel_sun   ! outer radius in solar radii
  real(kind=real64)    :: lnMaxT           ! running maximum of ln(T) over all zones
  real(kind=real64)    :: lnTCenter        ! ln(T) at the innermost zone (centre)
  real(kind=real64)    :: lnRhoCenter      ! ln(rho) at the innermost zone (centre)
  real(kind=real64)    :: model_mass       ! total mass from model metadata [M☉]
  real(kind=real64)    :: isoMass          ! isotope mass contribution for one zone [g]
  real(kind=real64) :: integratedMass ! numerically integrated total mass [g]
  real(kind=real64) :: thisMass       ! shell mass for current zone [g]
  real(kind=real64) :: thisIsoMass       ! isotopic shell mass for current zone [g]
  real(kind=real64) :: density        ! density for current zone [g/cm³]
  real(kind=real64) :: rinner, router ! inner/outer radius of current shell [cm]
  real(kind=real64) :: delradius      ! r_outer³ - r_inner³  [cm³]

  ! ---------------------------------------------------------------------------
  ! Isotope / element accounting
  ! ---------------------------------------------------------------------------
  integer :: num_isotopes              ! number of isotope columns in the model

  type(element_type) :: element
  ! ---------------------------------------------------------------------------
  ! Metadata lookup
  ! ---------------------------------------------------------------------------
  type(kv_pair) :: meta_mass_kv_pair   ! result of get_meta for "M/Msun"

  ! ---------------------------------------------------------------------------
  ! Output formatting buffers (tab-separated columns)
  ! ---------------------------------------------------------------------------
  character(len=:), allocatable :: reactOut   ! header row for network-comparison table
  character(len=:), allocatable :: dataOut    ! data row  for network-comparison table
  integer :: reactSize = 0
  integer :: dataSize = 0
  character(len=50)  :: tempData   ! temporary formatted field

  ! ---------------------------------------------------------------------------
  ! Logicals for output control
  ! ---------------------------------------------------------------------------
  
  logical :: do_isotopes_summary = .false.
  logical :: do_element_summary = .false.
  logical :: do_general_summary = .false.
  logical :: do_element_summary_to_si = .false.
  logical :: do_all = .false.
  
  ! ---------------------------------------------------------------------------
  ! List of model files to read and process
  ! ---------------------------------------------------------------------------
  
  character(len=256), allocatable :: model_list(:)
  integer :: num_models = 0
  
  character(len=8) :: isoID

  ! ==========================================================================
  ! 1. Parse command-line arguments
  ! ==========================================================================
  num_args = command_argument_count()
  allocate(args(num_args))

  do ix = 1, num_args
    call get_command_argument(ix, args(ix))
    if (args(ix)(1:1) == "-") then  ! look for run flags
      if (args(ix) == "-si") then
        do_element_summary_to_si = .true.
      end if
      if (args(ix) == "-i") then
        do_isotopes_summary = .true.
      end if
      if (args(ix) == "-e") then
        do_element_summary = .true.
      end if
      if (args(ix) == "-g") then
        do_general_summary = .true.
      end if
      if (args(ix) == "-all") then
        do_all = .true.
      end if
    else
      num_models = num_models + 1 ! count how many models there are to process
    end if
  end do
  
  allocate(model_list(num_models))
  
  ! read through the argument list and collect the paths to models to process. Hopefully the filename doesn't start with '-'
  num_models = 0
  do ix = 1, num_args
    if (args(ix)(1:1) .ne. "-") then
      num_models = num_models + 1
      model_list(num_models) = trim(args(ix))
    end if
  end do

  if (num_models == 0) then
    write(*,*) "Usage: mesa_mod_summary [-i] [-e] [-g] [-si] <model-file> <model-file> ..."
    write(*,*) "    -i: output integrated masses for each isotope"
    write(*,*) "    -si: output integrated masses for element up to si"
    write(*,*) "    -e: output integrated masses for each element"
    write(*,*) "    -g: output general summary"
    write(*,*) "    -all: output all isotopes or elements, regardless whether it is present"
    write(*,*) " if no options are selected, -g will be selected by default"
    stop
  end if
    ! ==========================================================================
    ! 2. Check and make sure at least one output option is selected. if none, do general summary
    ! ==========================================================================
  if (.not.do_element_summary_to_si.and..not.do_isotopes_summary.and..not.do_element_summary) then
    do_general_summary = .true.
  end if

    ! ==========================================================================
    ! 2. Initialize the isotope integrated mass table
    ! ==========================================================================
  call init_iso_mass()

  do nm = 1, num_models
  
    ! ==========================================================================
    ! 2. Reset the isotope integrated mass table
    ! ==========================================================================
    call reset_iso_mass() 
    ! ==========================================================================
    ! 2. Read model file
    ! ==========================================================================
    call read_mesa_model(model_list(nm), model)


    ! ==========================================================================
    ! 3. Extract total mass from model metadata
    !    The key "M/Msun" stores the mass as a quoted string; we do an internal
    !    read to obtain the numeric value.
    ! ==========================================================================
    call get_meta(model, "M/Msun", meta_mass_kv_pair)
    read(meta_mass_kv_pair%value, *) model_mass



    ! ==========================================================================
    ! 7. Extract scalar quantities from the model table
    !
    !    Row 1 is the outermost zone; row nrows is the innermost (centre).
    !    All spatial/thermodynamic values are stored as natural logarithms.
    ! ==========================================================================
    radius     = exp(model%table(1, 3))              ! outer radius [cm]
    radius_rel_earth = radius / earth_radius
    radius_rel_sun   = radius / solar_radius
    lnTCenter  = model%table(model%nrows, 2)         ! ln T at centre
    lnRhoCenter = model%table(model%nrows, 1)        ! ln rho at centre
    lnMaxT     = -1.0e9                              ! sentinel; will be updated


    ! ==========================================================================
    ! 8. Integrate mass zone by zone
    !
    !    Each zone is treated as a spherical shell.  The shell volume is
    !    computed from the outer and inner radii stored in adjacent rows:
    !
    !      V_shell = (4π/3) * (r_outer³ - r_inner³)
    !      M_shell = ρ * V_shell
    !
    !    The innermost zone has r_inner = 0 (stellar centre).
    !    The mass fraction for isotope k in zone irow is stored in
    !    table(irow, 6+k), so the isotope mass contribution is:
    !
    !      M_iso_shell = X_k * M_shell
    ! ==========================================================================
    integratedMass = 0.0d0

    do irow = 1, model%nrows

    ! Track the maximum temperature (in ln-space)
      if (model%table(irow, 2) > lnMaxT) then
        lnMaxT = model%table(irow, 2)
      end if

    ! Shell geometry
      density = exp(model%table(irow, 1))   ! [g/cm³]
      router  = exp(model%table(irow, 3))   ! outer edge of this shell [cm]
      if (irow == model%nrows) then
        rinner = 0.0d0                      ! centre of the star
      else
        rinner = exp(model%table(irow + 1, 3))  ! inner edge = next row's outer radius
      end if

    ! Shell mass  M = ρ * (4π/3) * (r_out³ - r_in³)
      delradius = router**3 - rinner**3
      thisMass  = density * (4.0d0 / 3.0d0) * PI * delradius

      integratedMass = integratedMass + thisMass

    ! Accumulate isotope and element masses using mass fractions from the table
      do ix = 1, model%num_isotopes
        isoMass = model%table(irow, 6 + ix) * thisMass
        isoID = model%col_header(6 + ix)
        if (isoID == "neut")  then 
          isoID = "n1"
        else
          if ((isoID == "prot").or.(isoID == "Prot"))  then 
            isoID = "H1"
          else
            isoID = title_case(model%col_header(6 + ix))
          end if
        end if
        call add_iso_mass(isoID,isoMass)
!        write (*,*) "Add iso mass to ", isoID, " of ", isoMass
      end do

    end do

    ! ==========================================================================
    ! 10. Print general stellar summary
    ! ==========================================================================
    if (do_general_summary) then
      write(*,*) "log Maximum temperature: ", log10(exp(lnMaxT)),  " log K"
      write(*,*) "log Central temperature: ", log10(exp(lnTCenter)), " log K"

    ! Display radius in solar or Earth radii depending on which is more natural
      if (radius_rel_earth > 9.95) then
        write(*,*) "Outer radius: ", radius, " cm (", radius_rel_sun, " R☉)"
      else
        write(*,*) "Outer radius: ", radius, " cm (", radius_rel_earth, " R⊕)"
      end if

      write(*,*) "Integrated Mass: ", integratedMass, " g (", &
             integratedMass / solar_mass, " M☉)"
      write(*,*) "Model mass (from metadata): ", model_mass, " M☉"

    ! Element-by-element breakdown
      do jx = 1, 85
        thisIsoMass = get_ele_mass_by_Z(jx)
        element = find_element_by_Z(jx)
        if (thisIsoMass > 0) then
          write(*,*) element%symbol, &
                 thisIsoMass / solar_mass, " M☉", &
                 thisIsoMass / integratedMass * 100.0, " %"
        end if
      end do


    ! ==========================================================================
    ! 11. Tab-separated output block for CO white-dwarf comparisons
    ! ==========================================================================
      write(*,*)
      write(*,*) "--- CO WD summary (tab-separated) ---"
      write(*,'(A, A, A, A, A)') "mass (solar)", achar(9), "radius (R⊕)", &
                               achar(9), "core density log(g/cm³)"
      write(*,'(F9.7, A, F13.7, A, F9.7)') &
        integratedMass / solar_mass, achar(9), &
        radius_rel_earth,            achar(9), &
        log10(exp(lnRhoCenter))
        
    end if! do general summary

    ! ==========================================================================
    ! 12. Tab-separated output block for nuclear-network comparisons
    !     Columns are added only for elements actually present in the network.
    ! ==========================================================================
    if (do_element_summary_to_si) then
      write(*,*)
      
      ! allocate reactOut
      reactSize = 12+25+19*6+10
      dataSize = 256+23+14*6+10
      allocate(character(len=reactSize) :: reactOut)
      allocate(character(len=dataSize) :: dataOut)

    ! Build header and data strings dynamically
      write(reactOut, '(A, A, A)') "mass (solar)", achar(9), "radius (R⊕)"
      write(dataOut,  '(F9.7, A, F13.7)') integratedMass / solar_mass, achar(9), &
                                        radius_rel_earth

      call printMassSolar(get_ele_mass_by_Z(2),tempData) ! He
      reactOut = trim(reactOut) // achar(9) // "He mass (solar)"
      dataOut  = trim(dataOut)  // achar(9) // trim(tempData)

      call printMassSolar(get_ele_mass_by_Z(6),tempData) ! C
      reactOut = trim(reactOut) // achar(9) // "C mass (solar)"
      dataOut  = trim(dataOut)  // achar(9) // trim(tempData)

      call printMassSolar(get_ele_mass_by_Z(8),tempData) ! O
      reactOut = trim(reactOut) // achar(9) // "O mass (solar)"
      dataOut  = trim(dataOut)  // achar(9) // trim(tempData)

      call printMassSolar(get_ele_mass_by_Z(10),tempData) ! Ne
      reactOut = trim(reactOut) // achar(9) // "Ne mass (solar)"
      dataOut  = trim(dataOut)  // achar(9) // trim(tempData)

      call printMassSolar(get_ele_mass_by_Z(12),tempData) ! Mg
      reactOut = trim(reactOut) // achar(9) // "Mg mass (solar)"
      dataOut  = trim(dataOut)  // achar(9) // trim(tempData)

      call printMassSolar(get_ele_mass_by_Z(14),tempData) ! Si
      reactOut = trim(reactOut) // achar(9) // "Si mass (solar)"
      dataOut  = trim(dataOut)  // achar(9) // trim(tempData)

      write(*,*) trim(reactOut)
      write(*,*) trim(dataOut)
      
      deallocate(reactOut)
      deallocate(dataOut)

    end if    

    ! ==========================================================================
    ! 13. Tab-separated output block for isotope integrated masses
    ! ==========================================================================
    if (do_isotopes_summary) then
    ! Build header and data strings dynamically
      reactSize = 12+25+19*n_isotopes+10
      dataSize = 256+23+14*n_isotopes+10
      allocate(character(len=reactSize) :: reactOut)
      allocate(character(len=dataSize) :: dataOut)

      if (num_models > 1) then
        write(reactOut, '(A, A, A, A, A)') "model path", achar(9), "mass (solar)", achar(9), "radius (R⊕)"
        write(dataOut,  '(A, A, F11.7, A, F13.7)') trim(model_list(nm)), achar(9), integratedMass / solar_mass, achar(9), &
                                      radius_rel_earth
      else 
        write(reactOut, '(A, A, A)') "mass (solar)", achar(9), "radius (R⊕)"
        write(dataOut,  '(F11.7, A, F13.7)') integratedMass / solar_mass, achar(9), &
                                      radius_rel_earth
      end if

      do ix = 1, n_isotopes
        isoMass = get_iso_mass_by_idx(ix)
        if ((isoMass > 0) .or. (do_all)) then
          call printMassSolar(isoMass,tempData)
          reactOut = trim(reactOut) // achar(9) // trim(isotope_table(ix)%id) // " mass (solar)"
          dataOut  = trim(dataOut)  // achar(9) // trim(tempData)
        end if
      end do

      if (nm == 1) then
        write(*,*) trim(reactOut)
      end if
      write(*,*) trim(dataOut)

      deallocate(reactOut)
      deallocate(dataOut)
    end if
 
    ! ==========================================================================
    ! 13. Tab-separated output block for isotope integrated masses
    ! ==========================================================================
    if (do_element_summary) then
    ! Build header and data strings dynamically
      reactSize = 12+25+16*n_elements+10
      dataSize = 256+23+14*n_elements+10
      allocate(character(len=reactSize) :: reactOut)
      allocate(character(len=dataSize) :: dataOut)

      if (num_models > 1) then
        write(reactOut, '(A, A, A, A, A)') "model path", achar(9), "mass (solar)", achar(9), "radius (R⊕)"
        write(dataOut,  '(A, A, F11.7, A, F13.7)') trim(model_list(nm)), achar(9), integratedMass / solar_mass, achar(9), &
                                        radius_rel_earth
      else 
        write(reactOut, '(A, A, A)') "mass (solar)", achar(9), "radius (R⊕)"
        write(dataOut,  '(F11.7, A, F13.7)') integratedMass / solar_mass, achar(9), &
                                      radius_rel_earth
      end if
      
      do ix = 1, n_elements
        isoMass = get_ele_mass_by_Z(ix)
        if ((isoMass > 0) .or. (do_all)) then
          call printMassSolar(isoMass,tempData)
          element = find_element_by_Z(ix)
          reactOut = trim(reactOut) // achar(9) // trim(element%symbol) // " mass (solar)"
          dataOut  = trim(dataOut)  // achar(9) // trim(tempData)
        end if
      end do

      if (nm == 1) then
        write(*,*) trim(reactOut)
      end if
      write(*,*) trim(dataOut)

      deallocate(reactOut)
      deallocate(dataOut)
    end if
  
    call destroy_mesa_model(model)
  end do
end program mesa_mod_summary
