!-----------------------------------------------------------------------
!Module: quantum
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! physics module
!!----------------------------------------------------------------------
!! Included subroutines:
!! sample box 
!! construct initial wave fnction 
!! evolve wave function 
!! construct time evolution matrix 
!! expectation values 
!! construct hamiltonian 
!! calculate analytic 
!! construct sho hamiltonian 
!! construct super matrix 1 
!! construct super matrix 2
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module quantum

use types
use linear_algebra, only : invert_matrix
implicit none
private
public sample_box, construct_initial_wavefunction, evolve_wave_function, construct_time_evolution_matrix, & 
expectation_values, construct_hamiltonian, calculate_analytic, construct_sho_hamiltonian

contains

!-----------------------------------------------------------------------
!! Subroutine: sample_box
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Fills an array with points along the potential from -L to L
!!----------------------------------------------------------------------
!! Input:
!! n_points     points from -L to L 
!! length       length of the potential 
!!----------------------------------------------------------------------
!! Output: 
!! dx           increment value  
!! x_vector     array containing points along the potential
!!----------------------------------------------------------------------
subroutine sample_box(x_vector, n_points, dx, length)
    implicit none 
    real(dp), allocatable, intent(out):: x_vector(:) 
    real(dp), intent(out) :: dx
    real(dp), intent(in) :: length
    integer, intent(in) :: n_points 
    integer :: i
    dx = 2*length/(n_points - 1)
    allocate(x_vector(1:n_points)) 
    !fill x array with points along the potential 
    !first element should be -L and last element should be L
    do i = 1, n_points 
        x_vector(i) = -length + dx * (i - 1) 
    end do 
    
end subroutine sample_box

!-----------------------------------------------------------------------
!! Subroutine: construct_initial_wavefunction
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Initializes wave function as Gaussian
!!----------------------------------------------------------------------
!! Input:
!! n_points     points along the potential 
!! x            x_vector array containg points along the potential 
!! width        initial width of the gaussian 
!! center       center of Gaussian 
!!----------------------------------------------------------------------
!! Output:
!! wave_function    array containing initial gaussian as a function of x 
!!----------------------------------------------------------------------
subroutine construct_initial_wavefunction(n_points, wave_function, x, width, center)
    implicit none
    real(dp), allocatable, intent(out) :: wave_function(:) 
    real(dp), intent(in) :: x(:), width, center
    real(dp) :: constant, exponential
    integer, intent(in) :: n_points 
    integer :: wave_points, i 
    wave_points = 2 * n_points 
    constant = (2 * pi * (width**2))**(-1.0_dp/4.0_dp)
    allocate(wave_function(1:wave_points)) 
    wave_function = 0.0_dp
    !Fill wave function array with value of gaussian at a specific x
    do i = 1, n_points 
        exponential = exp(-((x(i)-center)**2)/(4 * width**2))
        wave_function(i) = constant * exponential 
    end do 
    
end subroutine construct_initial_wavefunction 


!-----------------------------------------------------------------------
!! Subroutine: construct_hamiltonian
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Constructs Hamiltonian for potential of V(x) = 0
!!----------------------------------------------------------------------
!! Input:
!! dx   increment value 
!! n_points   number of points from -L to L 
!!----------------------------------------------------------------------
!! Output:
!! hamiltonian   Matrix containining kinetic and potential energy 
!!               operators
!!----------------------------------------------------------------------
subroutine construct_hamiltonian(dx, n_points, hamiltonian) 
real(dp), allocatable, intent(out) :: hamiltonian(:, :) 
real(dp), intent(in) :: dx
integer, intent(in) :: n_points 
integer :: i 

!Since V(x) = 0 Hamiltonian is just kinetic 
!energy terms. 
allocate(hamiltonian(1:n_points, 1:n_points)) 
hamiltonian = 0.0_dp
do i = 1, n_points 
    hamiltonian(i, i) = 1.0_dp/(dx**2) 
    hamiltonian(i, i + 1) = -0.5_dp/(dx**2) 
    hamiltonian(i + 1, i) = -0.5_dp/(dx**2)
end do 
end subroutine construct_hamiltonian



!-----------------------------------------------------------------------
!! Subroutine: construct_sho_hamiltonian
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Constructs Hamiltonian for harmonic oscillator
!!----------------------------------------------------------------------
!! Input:
!! dx   x increment value 
!! n_points   number of points from -L to L 
!! x_vector   array containing x value from -L to L 
!! k          wave number of packet for oscillator
!!----------------------------------------------------------------------
!! Output:
!! hamiltonian   matrix for hamiltonian operator of harmonic oscillator
!!----------------------------------------------------------------------
subroutine construct_sho_hamiltonian(dx, n_points, hamiltonian, x_vector, k) 
real(dp), allocatable, intent(out) :: hamiltonian(:, :) 
real(dp), intent(in) :: dx, x_vector(:), k 
real(dp), allocatable :: kinetic(:, :), potential(:, :)
integer, intent(in) :: n_points 
integer :: i 


allocate(hamiltonian(1:n_points, 1:n_points)) 
allocate(kinetic(1:n_points, 1:n_points)) 
allocate(potential(1:n_points, 1:n_points)) 
kinetic = 0.0_dp 
potential = 0.0_dp 
hamiltonian = 0.0_dp
!Kinetic energy doesn't change
do i = 1, n_points 
    kinetic(i, i) = 1.0_dp/(dx**2) 
    kinetic(i, i + 1) = -0.5_dp/(dx**2) 
    kinetic(i + 1, i) = -0.5_dp/(dx**2)
end do 
!Potential energy is 1/2 * k * x^2
do i = 1, n_points 
    potential(i, i) = 0.5_dp * k * (x_vector(i)**2) 
end do 
!Hamiltonian is defined as T + V 
hamiltonian = kinetic + potential 
end subroutine construct_sho_hamiltonian




!-----------------------------------------------------------------------
!! Subroutine: construct_time_evolution_matrix
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Constructs time evolution matrix
!!----------------------------------------------------------------------
!! Input:
!! hamiltonian      matrix containing hamiltonian operator 
!! n_points         points between -L and L 
!! delta_t          time increment value
!!----------------------------------------------------------------------
!! Output:
!! time_evolution   matrix that multiplies the wave function to evolve it 
!!                  forward in time
!!----------------------------------------------------------------------
subroutine construct_time_evolution_matrix(hamiltonian, time_evolution, n_points, delta_t)
    implicit none 
    real(dp), allocatable, intent(in) :: hamiltonian(:, :) 
    real(dp), allocatable, intent(out):: time_evolution(:, :) 
    real(dp), intent(in) :: delta_t
    real(dp), allocatable :: super_matrix_1(:, :), super_matrix_2_inv(:, :), super_matrix_2(:, :)
    integer, intent(in) :: n_points 
    integer :: time_points 
    
    time_points = 2 * n_points 
    !allocate arrays
    allocate(time_evolution(1:time_points, 1:time_points)) 
    allocate(super_matrix_1(1:time_points, 1:time_points)) 
    allocate(super_matrix_2_inv(1:time_points, 1:time_points)) 
    allocate(super_matrix_2(1:time_points, 1:time_points)) 
    !initialize arrays
    time_evolution = 0.0_dp 
    time_evolution = 0.0_dp
    super_matrix_1 = 0.0_dp
    super_matrix_2 = 0.0_dp 
    
    !call subroutines to construct matrices
    call construct_super_matrix_1(time_points, hamiltonian, delta_t, super_matrix_1) 
    print *, "constructed super matrix 1"
    call construct_super_matrix_2(time_points, hamiltonian, delta_t, super_matrix_2) 
    print *, "constructed super matrix 2" 
    !subroutine in linear algebra module to invert matrices
    call invert_matrix(super_matrix_2, super_matrix_2_inv) 
    print *, "inverted super matrix 2" 
    !use built in matrix multiplication
    time_evolution = matmul(super_matrix_2_inv, super_matrix_1)

end subroutine construct_time_evolution_matrix

!-----------------------------------------------------------------------
!! Subroutine: construct_super_matrix_1
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Constructs super matrix
!!----------------------------------------------------------------------
!! Input:
!! time_points    size of matrix 
!! hamiltonian    hamiltonian operator as a matrix 
!! delta_t        time increment value 
!!----------------------------------------------------------------------
!! Output:
!! super_matrix_1   super_matrix with identity and hamiltonian
!-----------------------------------------------------------------------

subroutine construct_super_matrix_1(time_points, hamiltonian, delta_t, super_matrix_1) 
implicit none 
real(dp), intent(out) :: super_matrix_1(:, :) 
real(dp), intent(in) :: delta_t, hamiltonian(:, :) 
integer, intent(in) :: time_points 
integer :: i, j, k 

 !Identity matrices
 do i = 1, time_points 
        super_matrix_1(i, i) = 1.0_dp 
 end do 

!Hamiltonian matrices multiplied by delta_t/2hbar
do i = 1, time_points/2 
        k = 1
        do j = (time_points/2) + 1, time_points  
            super_matrix_1(i, j) = hamiltonian(i, k) * (delta_t/(2)) 
            k = k + 1 
        end do 
    end do 

do j = 1, time_points/2 
    k = 1
        do i = (time_points/2) + 1, time_points 
            
            super_matrix_1(i, j) = hamiltonian(k, j) * (-delta_t/(2)) 
            
            k = k + 1 
            
        end do 
end do 

end subroutine construct_super_matrix_1


!-----------------------------------------------------------------------
!! Subroutine: construct_super_matrix_2
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Constructs second super matrix which will be inverted later
!!----------------------------------------------------------------------
!! Input: 
!! time_points    size of matrix 
!! hamiltonian    hamiltonian operator as a matrix 
!! delta_t        time increment value 
!!----------------------------------------------------------------------
!! Output:
!! super_matrix_2 super matrix containing identity and hamiltonian
!!----------------------------------------------------------------------

subroutine construct_super_matrix_2(time_points, hamiltonian, delta_t, super_matrix_2) 
implicit none 
real(dp), intent(out) :: super_matrix_2(:, :) 
real(dp), intent(in) :: delta_t, hamiltonian(:, :) 
integer, intent(in) :: time_points 
integer :: i, j, k 

 !Set up identity
 do i = 1, time_points 
        super_matrix_2(i, i) = 1.0_dp 
 end do 

!set up hamiltonian multiplied by delta_t/2hbar
do i = 1, time_points/2 
        k = 1
    do j = (time_points/2) + 1, time_points  
        super_matrix_2(i, j) = hamiltonian(i, k) * (-delta_t/(2)) 
        k = k + 1 
        end do 
 end do 

    do j = 1, time_points/2 
        k = 1
        do i = (time_points/2) + 1, time_points 
            super_matrix_2(i, j) = hamiltonian(k, j) *  (delta_t/(2)) 
            k = k + 1 
        end do 
    end do 

end subroutine construct_super_matrix_2


!-----------------------------------------------------------------------
!! Subroutine: evolve_wave_function
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Evolves wave function forward in time by multiplying it by the 
!! time evolution matrix
!!----------------------------------------------------------------------
!! Input:
!! time_evolution   time_evolution matrix 
!! wave_function    initial wave_function 
!! n_points         points from -L to L 
!!----------------------------------------------------------------------
!! Output:
!! time_wave_function   wave function at every time step as a matrix 
!! density              matrix with probability densities at every time steo
!!----------------------------------------------------------------------
subroutine evolve_wave_function(time_evolution, time_wave_function, wave_function, n_points, n_steps, density)
    implicit none
    real(dp), intent(in) :: time_evolution(:, :)
    real(dp), intent(inout) :: wave_function(:)
    real(dp), allocatable, intent(out) :: time_wave_function(:, :), density(:, :) 
    integer, intent(in) :: n_points, n_steps 
    integer :: i, j, k 
    allocate(time_wave_function(1:n_points, 1:n_steps + 1)) 
    allocate(density(1:n_points, 1:n_steps + 1)) 
    !Initialize first columns to be initial wave function values
    do i = 1, n_points
        time_wave_function(i, 1) = wave_function(i) 
        density(i, 1) = wave_function(i) **2 + wave_function(i + n_points)**2
    end do 
    !evolve forward in time my multiplying n_steps + 1 times 
    do i = 2, n_steps + 1 
       wave_function = matmul(time_evolution, wave_function)
       do j = 1, n_points 
       density(j, i) = wave_function(j)**2 + wave_function(j + n_points)**2
       time_wave_function(j, i) = wave_function(j) 
       end do 
    end do  
end subroutine evolve_wave_function

!-----------------------------------------------------------------------
!! Subroutine: calculate_norm
!-----------------------------------------------------------------------
!! By Nathan Crawford
!!
!! Calculates normalization constant at each time step
!!----------------------------------------------------------------------
!! Input:
!! density  array containing probability densities 
!! dx       x increment value 
!! n_steps  time steps 
!! n_points points from -L to L
!!----------------------------------------------------------------------
!! Output:
!! norm     array containing normalization constant at each time step
!!----------------------------------------------------------------------
subroutine calculate_norm(density, dx, norm, n_steps, n_points)  
implicit none 
real(dp), intent(out) :: norm(:) 
real(dp), intent(in) :: dx, density(:, :) 
real(dp) :: summation
integer, intent(in) :: n_steps, n_points 
integer :: i, j, k

!we can compute normalization by summing 
!the probability densities at a time and multiplying 
!by dx
do i = 1, n_steps + 1 
    summation = 0.0_dp 
    
        do j = 1, n_points 
            
            summation = summation + (density(j, i)  * dx)  
        end do 
    
    norm(i) = summation 
    
end do  
end subroutine calculate_norm

!-----------------------------------------------------------------------
!! Subroutine: expectation_values
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculates expectation values of position and standard deviation
!!----------------------------------------------------------------------
!! Input:
!! norm     array containing normalization constants 
!! n_steps  time steps 
!! n_points points from -L to L 
!! density  matrix containing probability densities at each time step 
!! dx       x increment value 
!! x_vector array with x values from -L to L
!!----------------------------------------------------------------------
!! Output:
!! position         array with expectation value of position at each time 
!! position_squared array with expectation value of x^2 at each time 
!! sigma            array with standard deviation at each time
!!----------------------------------------------------------------------
subroutine expectation_values(norm, position, sigma, position_squared, n_steps, n_points, density, dx, x_vector)
    implicit none
    real(dp), allocatable, intent(out) :: position(:), sigma(:), position_squared(:), norm(:) 
    real(dp), intent(in) :: dx, density(:, :), x_vector(:) 
    real(dp) :: summation 
    integer, intent(in) :: n_steps, n_points 
    integer :: i, j, k 
    !allocate arrays
    allocate(position(1:n_steps + 1)) 
    allocate(sigma(1:n_steps + 1)) 
    allocate(position_squared(1:n_steps + 1)) 
    allocate(norm(1:n_steps + 1)) 
    !call subroutine for normalization constants
    call calculate_norm(density, dx, norm, n_steps, n_points) 
    !calculate expectation value of position 
    do i = 1, n_steps + 1 
        summation = 0.0_dp
        do j = 1, n_points 
            summation = summation + (x_vector(j) * density(j, i) * dx) 
            
        end do 
        position(i) = summation/norm(i)
    end do
     !calculate expectation value of x^2
     do i = 1, n_steps + 1 
        summation = 0.0_dp
        do j = 1, n_points 
            summation = summation + (x_vector(j)**2 * density(j, i) * dx) 
        end do 
        position_squared(i) = summation/norm(i)
    end do
    !calculate standard deviation
    do i = 1, n_steps + 1 
        sigma(i) = sqrt(position_squared(i) - (position(i)**2))
    end do 
end subroutine expectation_values 

!-----------------------------------------------------------------------
!! Subroutine: calculate_analytic
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! Calculates analytic expressions for sigma and position
!!----------------------------------------------------------------------
!! Input:
!! k        wave packet number 
!! center   center of gaussian 
!! width    width of gaussian 
!! delta_t  time increment value 
!! n_steps  number of time steps 
!!----------------------------------------------------------------------
!! Output:
!! analytic_position    values of analytic expression for position 
!! analytic_sigma       values of analytic expression for standard deviation
!!----------------------------------------------------------------------
subroutine calculate_analytic(analytic_position, analytic_sigma, k, center, width, delta_t, n_steps) 
implicit none 
real(dp), intent(out), allocatable :: analytic_position(:), analytic_sigma(:) 
real(dp), intent(in) :: k, center, width, delta_t 
real(dp) :: t 
integer, intent(in) :: n_steps 
integer :: i 
!allocate arrays
allocate(analytic_position(1:n_steps + 1)) 
allocate(analytic_sigma(1:n_steps + 1)) 
t = 0.0_dp 
!compute analytic expressions for standard deviation and 
!position of harmonic oscillator
do i = 1, n_steps + 1 
    analytic_sigma(i) = sqrt(width**2 + (t**2/(4*width**2))) 
    analytic_position(i) = center * cos(sqrt(k)*t) 
    t = t + delta_t
end do
end subroutine calculate_analytic 
  
end module quantum
