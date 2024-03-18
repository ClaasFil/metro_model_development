import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
# Assuming the numerical data provided in the question is stored in a file, let's read the file and parse the data into arrays.
# The file is expected to have three lines:
# 1. Numerical 2nd derivative values
# 2. Analytically calculated 2nd derivative values
# 3. Error values


# Let's read the file and parse the data.

file_name = 'task3_res.txt'  # Change this to the actual file name
file_path = os.path.join(os.getcwd(), file_name)

# Initialize empty lists to store the data
numerical_derivatives = []
analytical_derivatives = []
errors = []

with open(file_path, 'r') as file:
    lines = file.readlines()
    
    # Process each line
    count = 0
    for line in lines:
        # Split the line into words
        parts = line.split()
        
        # Skip the first four words ("Numerical 2nd derivative:", etc.) and convert the rest to floats
        if count == 1:
            numerical_derivatives = [float(value) for value in parts]
        elif count == 3:
            analytical_derivatives = [float(value) for value in parts]
        elif count == 5:
            errors = [float(value) for value in parts]

        count += 1

numerical_derivatives = np.array(numerical_derivatives)
analytical_derivatives = np.array(analytical_derivatives)
errors = np.array(errors)

# Index array for x-axis
index_array = np.arange(len(numerical_derivatives))

#interpolate the analytical derivatives
x = np.arange(0, len(analytical_derivatives))
y = analytical_derivatives
x_new = np.arange(0, len(analytical_derivatives), 0.1)

# Perform cubic interpolation on the analytical derivatives
cs = CubicSpline(x, y)
y_new = cs(x_new)


print("Numerical 2nd derivative values:", numerical_derivatives)
print("Analytically calculated 2nd derivative values:", analytical_derivatives)
print("Error values:", errors)

# Plot numerical and analytical derivatives
plt.figure(figsize=(14, 6))

plt.subplot(1, 2, 1)  # 1 row, 2 columns, 1st subplot

plt.scatter(index_array, numerical_derivatives, label='Numerical Derivatives', marker='x')
plt.scatter(index_array, analytical_derivatives, label='Analytical Derivatives', marker='x')
plt.plot(x_new, y_new, label='Interpolated Analytical Derivatives', linestyle='dashed')
plt.title('Numerical vs Analytical Derivatives')
plt.xlabel('Index')
plt.ylabel('Value')
plt.legend()

# Plot the errors
plt.subplot(1, 2, 2)  # 1 row, 2 columns, 2nd subplot
plt.plot(errors[1:-1], label='Error', color='red', marker='o')
plt.title('Error in Derivatives')
plt.xlabel('Datapoint Index')
plt.ylabel('Error Value')
plt.legend()

plt.tight_layout()
plt.show()

