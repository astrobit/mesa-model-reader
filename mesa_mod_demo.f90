
! =============================================================================
! Example driver program — compile together with the module above, e.g.:
!   gfortran -o demo mesa_mod_reader.f90
! =============================================================================
program mesa_mod_demo
  use mesa_mod_reader
  implicit none

  integer :: num_args
  character(len=256), dimension(:), allocatable :: args
  type(mesa_model) :: m
  integer :: i, icol

  num_args = command_argument_count()
  allocate(args(num_args))  ! I've omitted checking the return status of the allocation 

  do i = 1, num_args
      call get_command_argument(i,args(i))
         ! now parse the argument as you wish
  end do

  if (num_args == 0) then
    write(*,*) "Specify file to read."
    STOP
  end if

  call read_mesa_model(args(1), m)

  if (m%ncols < 0) then
    write(*,'(a)') 'Failed to read model.'
    stop 1
  end if

  write(*,'(a,i0)') 'ncols  : ', m%ncols
  write(*,'(a,i0)') 'nrows  : ', m%nrows

  write(*,'(/,a)') '--- Pre-table metadata ---'
  do i = 1, m%n_meta
    write(*,'(2a,a,a)') '  ', trim(adjustl(m%meta(i)%label)), ' = ', &
                          trim(m%meta(i)%value)
!    if (len_trim(m%meta(i)%label2) > 0) &
!      write(*,'(2a,a,a)') '  ', trim(adjustl(m%meta(i)%label2)), ' = ', &
!                            trim(m%meta(i)%value2)
  end do

  write(*,'(/,a)') '--- Column headers ---'
  do icol = 1, m%ncols
    write(*,'(i4,2a)') icol, ': ', trim(m%col_header(icol))
  end do

  write(*,'(/,a)') '--- First 3 data rows ---'
  do i = 1, min(3, m%nrows)
    write(*,'(a,i0)') '  row_number = ', m%row_number(i)
    do icol = 1, m%ncols
      write(*,'(4x,a,a,1pe25.16)') trim(m%col_header(icol)), ' = ', &
                                    m%table(i, icol)
    end do
  end do

  write(*,'(/,a)') '--- Post-table metadata ---'
  do i = 1, m%n_post_meta
    if (len_trim(m%post_meta(i)%label) > 0) &
      write(*,'(2a,a,a)') '  ', trim(adjustl(m%post_meta(i)%label)), ' = ', &
                            trim(m%post_meta(i)%value)
  end do

  call destroy_mesa_model(m)

end program mesa_mod_demo
