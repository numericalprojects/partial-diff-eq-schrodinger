! Program: schrodinger
! By: Nathan Crawford
!-----------------------------------------------------------------------------
! This program solves the Time Dependent Schrodinger Equation with the 
! Crank Nicolson method by first taking input from a namelist, initializing the 
! wave function, constructing the time evolution matrix and evolving the wave 
! function forward in time. 
! It then writes the resulting probability densities and expectation values 
! into 2 separate files.
!-----------------------------------------------------------------------------
program schrodinger 

use types
use read_write, only : read_input, write_time_evolution, write_expectation_values, & 
                        write_time_evolution_sho, write_expectation_values_sho
use quantum, only : sample_box, construct_initial_wavefunction, construct_time_evolution_matrix, &
    evolve_wave_function, expectation_values, construct_hamiltonian, calculate_analytic, & 
    construct_sho_hamiltonian

implicit none

real(dp) :: length, delta_t, width, center, k_oscillator, dx
integer :: n_points, n_steps
character(len=1024) :: time_file, density_file 
real(dp), allocatable :: x_vector(:), density(:, :) !will be of size n_points.
real(dp), allocatable :: wave_function(:)! will be of size 2*n_points.
real(dp), allocatable :: evolution_matrix(:,:), hamiltonian(:, :) !will be of size 2*n_points by 2*n_points
real(dp), allocatable :: time_wave_function(:,:) !will be of size n_points by n_steps + 1 (the +1 is so that you can store the t=0 value)
real(dp), allocatable :: norm(:), position(:), position_squared(:), sigma(:), analytic_sigma(:), analytic_position(:) !all of size n_steps + 1 

call read_input(length, n_points, n_steps, delta_t, width, center, k_oscillator, time_file, density_file) 
print *, "read input"
call sample_box(x_vector, n_points, dx, length) 
call construct_initial_wavefunction(n_points, wave_function, x_vector, width, center) 
print *, "constructed wave"
call construct_hamiltonian(dx, n_points, hamiltonian) 
print *, "made hamil"
call construct_time_evolution_matrix(hamiltonian, evolution_matrix, n_points, delta_t) 
print *, "constructed time"
call evolve_wave_function(evolution_matrix, time_wave_function, wave_function, n_points, n_steps, density) 
print *, "evolved wave function"
call expectation_values(norm, position, sigma, position_squared, n_steps, n_points, density, dx, x_vector) 
print *, "calculated expectation values"
call calculate_analytic(analytic_position, analytic_sigma, k_oscillator, center, width, delta_t, n_steps) 
print *, "calculated analytic functions"
call write_time_evolution(density, x_vector, n_points, n_steps, density_file, wave_function)
call write_expectation_values(position, position_squared, analytic_sigma, sigma, time_file, n_steps, delta_t) 



call construct_initial_wavefunction(n_points, wave_function, x_vector, width, center) 
call construct_sho_hamiltonian(dx, n_points, hamiltonian, x_vector, k_oscillator) 
call construct_time_evolution_matrix(hamiltonian, evolution_matrix, n_points, delta_t) 
call evolve_wave_function(evolution_matrix, time_wave_function, wave_function, n_points, n_steps, density) 
call expectation_values(norm, position, sigma, position_squared, n_steps, n_points, density, dx, x_vector) 
call calculate_analytic(analytic_position, analytic_sigma, k_oscillator, center, width, delta_t, n_steps) 
call write_time_evolution_sho(density, x_vector, n_points, n_steps, density_file, wave_function)
call write_expectation_values_sho(position, position_squared, sigma, time_file, n_steps, delta_t, analytic_sigma, analytic_position) 
end program schrodinger
