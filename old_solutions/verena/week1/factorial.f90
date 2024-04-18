
program FactorialCalculator
	
	implicit none
	integer :: n, i
	integer :: factorial

	! We want to ask for  positive integer
	print *, "Please enter a positive integer:"
	read(*, *) n

	! We want to check if the integer is positive
	if (n <= 0) then
		print *, "Error: Number must be a positive integer."
	else
		! Calculate factorial
		factorial = 1
		do i = 1, n
			factorial = factorial * i
		end do

		! Print the result
		print *, "Factorial of ", n, " is ", factorial
	end if
end program FactorialCalculator

