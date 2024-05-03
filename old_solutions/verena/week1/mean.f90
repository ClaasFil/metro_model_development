program MeanStandardDeviation

	implicit none
	integer :: n, i
	real :: number, sum, sum_of_squares, mean, variance, std_deviation

	sum = 0.0
	sum_of_squares = 0.0

	! Ask for how many numbers
	print *, "How many numbers do you want to enter?"
	read (*, *) n

	! Read the numbers one by one and calculate sum of squares
	print *, "Enter the numbers one by one and hit enter after each number."
	do i = 1, n
		read(*, *) number
		sum = sum + number
		sum_of_squares = sum_of_squares + number**2
	end do

	! Calculate mean, variance and standard deviation
	mean = sum / real(n)
	variance = sum_of_squares / real(n) - mean**2
	std_deviation = sqrt(variance)

	! Output results
	print *, "Mean: ", mean
	print *, "Standard Deviation: ", std_deviation

end program MeanStandardDeviation

