program MeanArray

	implicit none
	integer :: n, i
	real, allocatable :: numbers(:)
	real :: sum_numbers, sum_of_squares, mean, variance, std_deviation
	
	print *, "How many numbers do you want to enter?"
	read(*, *) n
	
	! Check if number is negative
	do while (n<=0)
		print *, "Error: Please enter a positive number!"
		read(*, *) n
	end do

	allocate(numbers(n))

	print *, "Enter the numbers one by one: "
	read(*, *) numbers

	sum_numbers = sum(numbers)
	sum_of_squares = sum(numbers**2)
	mean = sum_numbers / real(n)
	variance = sum_of_squares / real(n) - mean**2
	std_deviation = sqrt(variance)

	print *, "Mean: ", mean
	print *, "Standard Deviation: ", std_deviation

end program MeanArray

