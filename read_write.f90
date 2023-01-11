!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! These subroutines take input from the user and write results into files
!!----------------------------------------------------------------------
!! Included subroutines:
!! read_input
!! write_time_evolution 
!! write_expectation_values 
!! write_time_evolution_sho 
!! write_expectation_values_sho
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module read_write
use types

implicit none

private
public :: read_input, write_time_evolution, write_expectation_values, & 
          write_time_evolution_sho, write_expectation_values_sho


contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Gets input from user from namelist, if no namelist then uses the 
!! default
!!----------------------------------------------------------------------
!! Input:
!! 
!!----------------------------------------------------------------------
!! Output:
!! length           length of the box 
!! n_points         points from -L to L 
!! n_steps          time steps 
!! delta_t          time increment value 
!! width            width of gaussian 
!! center           center of gaussian 
!! k_oscillator     wave number of packet 
!! time_file        file for expectation values 
!! density_file     file for probability densities
!!----------------------------------------------------------------------
subroutine read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator &
    , time_file, density_file)
    implicit none
    real(dp), intent(out) :: length, delta_t, width, center, k_oscillator
    integer, intent(out) :: n_points, n_steps
    character(len=*) :: time_file, density_file 
    integer :: ierror, file_unit, n_arguments 
    character(len = 200) :: namelist_file
    logical:: file_exists
    
    
    
    !namelist setup
    namelist /integration/ length, n_points, n_steps, delta_t
    namelist /wave_function/ width, center
    namelist /oscillator/ k_oscillator
    namelist /output/ time_file, density_file

    length = 5._dp
    n_points = 100
    n_steps = 100
    delta_t = 0.05_dp
    width = 0.5_dp
    center = 2.0_dp
    k_oscillator = 0.0_dp
    time_file = 'time_results.dat'
    density_file = 'density_results.dat' 
    
    n_arguments = command_argument_count()
    !namelist error trap
    if (n_arguments == 1) then
        call get_command_argument(1, namelist_file)
        inquire(file = trim(namelist_file), exist = file_exists)
        if (file_exists) then
            open(newunit=file_unit, file=namelist_file)
            read(file_unit, nml=integration, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading box namelist"
                stop
            endif
            read(file_unit, nml=wave_function, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading charge_distribution namelist"
                stop
            endif
            read(file_unit, nml=oscillator, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading sampling namelist"
                stop
            endif
            read(file_unit, nml=output, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading output namelist"
                stop
            endif
        else
            print*, namelist_file, 'not found'
            stop
        endif
    elseif (n_arguments /= 0) then
        print*, 'Incorrect number of arguments. Program takes either 0 or 1 argument only'
        print*, 'See details in README.md'
        stop
    endif
    print *, "greetings user, this program solves the time dependent" 
    print *, "schrodinger equation for 0 potential and harmonic oscillator."

end subroutine read_input


!-----------------------------------------------------------------------
!! Subroutine: write_time_evolution
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! writes densities as a function of x into a file
!!----------------------------------------------------------------------
!! Input:
!! density          probability density matrix 
!! x_vector         sampled points along the box 
!! n_points         points between -L and L 
!! n_steps          number of time steps 
!! density_file     file for output results 
!!  
!!----------------------------------------------------------------------
!! Output:
!!
!!----------------------------------------------------------------------
subroutine write_time_evolution(density, x_vector, n_points, n_steps, density_file, wave_function)
    implicit none
    real(dp), intent(in) :: wave_function(:), x_vector(:) 
    real(dp), allocatable :: density(:, :)
    character(len=1024), intent(in) :: density_file
    integer, intent(in) :: n_points, n_steps 
    integer :: i, j, k, file_unit  
    open(unit = file_unit, file = density_file)
    !we don't want all the values because that is a large file so we are only writing 
    !snapshots
    write(file_unit, *) "x ", "density ", "density ", "density ", "density "
    do i = 1, n_points 
        write(file_unit, *) x_vector(i), density(i, 1), density(i, (n_steps + 1)/3), & 
        density(i, 2*(n_steps + 1)/3), density(i, n_steps + 1)
    end do  
    close(file_unit)
    
end subroutine write_time_evolution

!-----------------------------------------------------------------------
!! Subroutine: write_expectation_values
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! writes expectation values as a function of time
!!----------------------------------------------------------------------
!! Input:
!! position         expectation value of position array 
!! position_squared expectation value of position_squared array 
!! analytic_sigma   array with sigma calculated analytically 
!! sigma            array with sigma calculated numerically 
!! time_file        file for output 
!! n_steps          time steps 
!! delta_t          time increment value
!!----------------------------------------------------------------------
!! Output:
!! 
!!----------------------------------------------------------------------
subroutine write_expectation_values(position, position_squared, analytic_sigma, sigma, time_file, n_steps, delta_t)
    implicit none
    real(dp), intent(in) :: position(:), position_squared(:), sigma(:), delta_t 
    real(dp), intent(in) :: analytic_sigma(:)
    real(dp) :: t
    character(len=1024), intent(in) :: time_file 
    integer, intent(in) :: n_steps 
    integer :: i, j, k, file_unit 
    open(unit = file_unit, file = time_file) 
    write(file_unit, *) "time ", "position ", "sigma ", "analytic sigma"
    t = 0.0_dp
    !expecation values as a function of time
    do i = 1, n_steps + 1 
        write(file_unit, *) t, position(i), position_squared(i), sigma(i), analytic_sigma(i)
        t = t + delta_t
    end do
    close(file_unit)
end subroutine write_expectation_values 


!-----------------------------------------------------------------------
!! Subroutine: write_time_evolution_sho
!-----------------------------------------------------------------------
!! By:
!!
!! Give an explanation of what the subroutine does
!!----------------------------------------------------------------------
!! Input:
!! density          probability density matrix 
!! x_vector         sampled points along the box 
!! n_points         points between -L and L 
!! n_steps          number of time steps 
!! density_file     file for output results 
!!----------------------------------------------------------------------
!! Output:
!!
!!----------------------------------------------------------------------
subroutine write_time_evolution_sho(density, x_vector, n_points, n_steps, density_file, wave_function)
    implicit none
    real(dp), intent(in) :: wave_function(:), x_vector(:) 
    real(dp), allocatable :: density(:, :)
    character(len=1024), intent(inout) :: density_file
    integer, intent(in) :: n_points, n_steps 
    integer :: i, j, k, file_unit  
    density_file = 'density_results_sho.dat' 
    open(unit = file_unit, file = density_file)
    
    write(file_unit, *) "x ", "density ", "density ", "density ", "density "
    do i = 1, n_points 
        write(file_unit, *) x_vector(i), density(i, 1), density(i, (n_steps + 1)/3), & 
        density(i, 2*(n_steps + 1)/3), density(i, n_steps + 1)
    end do  
    close(file_unit)
    
end subroutine write_time_evolution_sho

!-----------------------------------------------------------------------
!! Subroutine: write_expectation_values_sho
!-----------------------------------------------------------------------
!! By:
!!
!! Give an explanation of what the subroutine does
!!----------------------------------------------------------------------
!! Input:
!! position         expectation value of position array 
!! position_squared expectation value of position_squared array 
!! analytic_sigma   array with sigma calculated analytically 
!! sigma            array with sigma calculated numerically 
!! time_file        file for output 
!! n_steps          time steps 
!! delta_t          time increment value 
!! analytic_position array containing values for position calculated analytically
!!----------------------------------------------------------------------
!! Output:
!!
!!----------------------------------------------------------------------
subroutine write_expectation_values_sho(position, position_squared, sigma, time_file, n_steps, delta_t, & 
analytic_sigma, analytic_position)
    implicit none
    real(dp), intent(in) :: position(:), position_squared(:), sigma(:), delta_t, analytic_sigma(:) 
    real(dp), intent(in) :: analytic_position(:)
    real(dp) :: t
    character(len=1024), intent(inout) :: time_file 
    integer, intent(in) :: n_steps 
    integer :: i, j, k, file_unit 
    time_file = 'time_results_sho.dat' 
    t = 0.0_dp
    open(unit = file_unit, file = time_file) 
    write(file_unit, *) "time ", "position ", "sigma ", "analytic sigma ", "analytic position "
    do i = 1, n_steps + 1    
        write(file_unit, *) t, position(i), position_squared(i), sigma(i), analytic_sigma(i), analytic_position(i)
        t = t + delta_t
    end do
    close(file_unit) 
end subroutine write_expectation_values_sho

end module read_write
