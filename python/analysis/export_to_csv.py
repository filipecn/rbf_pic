import numpy as np
import csv


x_values = np.load("x_values.numpy.npy")
y_cond = np.load("y_cond.numpy.npy")
y_iterations = np.load("y_iterations.numpy.npy")
matrix_size = np.load("matrix_size.numpy.npy")

filename = 'regular_simulation.csv'
with open(filename, 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(["frame", "condition number", "number of iterations", "matrix size"])
    for item in zip(x_values, y_cond, y_iterations, matrix_size):
        writer.writerow(item)


# load the results - adaptive
x_values_adaptive = np.load("x_values.adaptive.npy")
y_cond_adaptive = np.load("y_cond.adaptive.npy")
y_iterations_adaptive = np.load("y_iterations.adaptive.npy")
matrix_size_adaptive = np.load("matrix_size.adaptive.npy")


filename = 'adaptive_simulation.csv'
with open(filename, 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(["frame", "condition number", "number of iterations", "matrix size"])
    for item in zip(x_values_adaptive, y_cond_adaptive, y_iterations_adaptive, matrix_size_adaptive):
        writer.writerow(item)

