import plotly
plotly.tools.set_credentials_file(username='jehutymax', api_key='KM31JtOX1x24IoA5EQzr')
import plotly.plotly as py
import plotly.graph_objs as go
import plotly.io as pio

import numpy as np
import scipy.io
import scipy.sparse.linalg
import math
import glob

# Load Matrix Market files generated during furoo simulations at the pressure solve stage.


# def compute_sparse_condition_number(matrix):
#     norm = scipy.sparse.linalg.norm(matrix)
#     norm_inverse = scipy.sparse.linalg.norm(scipy.sparse.linalg.inv(matrix))
#     return norm * norm_inverse
#
#
# def load_matrices(pattern):
#     path = "/home/rafael/Projects/bitbucket/furoo/build-release/tests/simulations"
#     buffer = glob.glob(path + "/" + pattern)
#     # remove rhs vectors
#     files = [s for s in buffer if "rhs" not in s]
#     entries = {}  # entries should be a list [num_iterations, condition_number, filename, matrix_size]
#
#     for file in files:
#         num_frame = int(file.split("_")[1])
#         num_iterations = int(file.split("_")[2].split("=")[1].split(".")[0])
#         matrix = scipy.io.mmread(file).tocsc()
#         cond_number = compute_sparse_condition_number(matrix)
#         num_rows = matrix.shape[0]
#         result = [num_iterations, cond_number, file, num_rows]
#         entries[num_frame] = result
#     return entries
#
#
# # load all files for regular simulation
# regular_entries = load_matrices("regular*.mtx")
# # adaptive_entries = load_matrices("adaptive*.mtx")
# #
# # # print(len(regular_files), regular_files)
# # # print(regular_entries)
# # print(adaptive_entries)
# x_values = np.zeros(len(regular_entries))
# y_cond = np.zeros(len(regular_entries))
# y_iterations = np.zeros(len(regular_entries))
# matrix_size = np.zeros(len(regular_entries))
# #
# # x_values = np.zeros(len(adaptive_entries))
# # y_cond = np.zeros(len(adaptive_entries))
# # y_iterations = np.zeros(len(adaptive_entries))
# # matrix_size = np.zeros(len(adaptive_entries))
# #
# counter = 0
# for frame in range(10, 300, 10): #regular_entries:
#     x_values[counter] = frame
#     y_cond[counter] = regular_entries[frame][1]
#     y_iterations[counter] = regular_entries[frame][0]
#     matrix_size[counter] = regular_entries[frame][3]
#     counter = counter + 1
# #
# # counter = 0
# # for frame in range(10, 300, 10): #adaptive_entries:
# #     x_values[counter] = frame
# #     y_cond[counter] = adaptive_entries[frame][1]
# #     y_iterations[counter] = adaptive_entries[frame][0]
# #     matrix_size[counter] = adaptive_entries[frame][3]
# #     counter = counter + 1
#
# # # save the results
# np.save("x_values.numpy", x_values)
# np.save("y_cond.numpy", y_cond)
# np.save("y_iterations.numpy", y_iterations)
# np.save("matrix_size.numpy", matrix_size)
# # save the results - adaptive
# # np.save("x_values.adaptive", x_values)
# # np.save("y_cond.adaptive", y_cond)
# # np.save("y_iterations.adaptive", y_iterations)
# # np.save("matrix_size.adaptive", matrix_size)

# # load the results
x_values = np.load("x_values.numpy.npy")
y_cond = np.load("y_cond.numpy.npy")
y_iterations = np.load("y_iterations.numpy.npy")
matrix_size = np.load("matrix_size.numpy.npy")
# load the results - adaptive
x_values_adaptive = np.load("x_values.adaptive.npy")
y_cond_adaptive = np.load("y_cond.adaptive.npy")
y_iterations_adaptive = np.load("y_iterations.adaptive.npy")
matrix_size_adaptive = np.load("matrix_size.adaptive.npy")

# trace1 = go.Scatter(
#     x=x_values,
#     y=y_cond,
#     mode="lines+markers",
#     name="Condition Number"
# )
#
# trace2 = go.Scatter(
#     x=x_values,
#     y=y_iterations,
#     mode="lines+markers",
#     name="Number of Iterations",
#     yaxis='y2'
# )
#
# trace3 = go.Scatter(
#     x=x_values,
#     y=matrix_size,
#     mode="lines+markers",
#     name="Matrix Size",
#     yaxis='y3'
# )
#
# data = [trace1, trace2, trace3]
# layout = go.Layout(
#     title="Pressure Matrix Analysis for 2D Regular Grid Simulation",
#     # title="Pressure Matrix Analysis for 2D Adaptive Grid Simulation",
#     yaxis=dict(
#         title="Condition Number"
#     ),
#     yaxis2=dict(
#         title="Number of Iterations",
#         titlefont=dict(
#             color='rgb(148, 103, 189)'
#         ),
#         tickfont=dict(
#             color='rgb(148, 103, 189)'
#         ),
#         anchor='x',
#         overlaying='y',
#         side='right',
#         position=0.85
#     ),
#     yaxis3=dict(
#             title="Matrix Size",
#             titlefont=dict(
#                 color='rgb(189, 103, 148)'
#             ),
#             tickfont=dict(
#                 color='rgb(189, 103, 148)'
#             ),
#             anchor='free',
#             overlaying='y',
#             side='right'
#     ),
#     legend=dict(orientation="h")
# )
#
# fig = go.Figure(data=data, layout=layout)
# py.plot(fig)

# trace1 = go.Scatter(
#     x=x_values,
#     y=y_cond,
#     mode="lines",
#     name="Regular Grid"
# )
#
# trace2 = go.Scatter(
#     x=x_values,
#     y=y_cond_adaptive,
#     mode="lines",
#     name="Adaptive Grid"
# )
#
#
# data = [trace1, trace2]
# layout = go.Layout(
#     title="PPE Matrix Condition Number",
#     yaxis=dict(
#         title="Condition Number"
#     ),
#     xaxis=dict(
#         title="Frame"
#     ),
#     # legend=dict(orientation="h")
# )
#
# fig = go.Figure(data=data, layout=layout)
# pio.write_image(fig, "ppe_cond_number.pdf")
# py.plot(fig)

# trace1 = go.Scatter(
#     x=x_values,
#     y=y_iterations,
#     mode="lines",
#     name="Regular Grid"
# )
#
# trace2 = go.Scatter(
#     x=x_values,
#     y=y_iterations_adaptive,
#     mode="lines",
#     name="Adaptive Grid"
# )
#
#
# data = [trace1, trace2]
# layout = go.Layout(
#     title="Number of Iterations Required to Solve PPE System",
#     yaxis=dict(
#         title="Number of Iterations"
#     ),
#     xaxis=dict(
#         title="Frame"
#     ),
#     # legend=dict(orientation="h")
# )
#
# fig = go.Figure(data=data, layout=layout)
# pio.write_image(fig, "ppe_iterations.pdf")
# py.plot(fig)

trace1 = go.Scatter(
    x=x_values,
    y=matrix_size,
    mode="lines",
    name="Regular Grid"
)

trace2 = go.Scatter(
    x=x_values,
    y=matrix_size_adaptive,
    mode="lines",
    name="Adaptive Grid"
)


data = [trace1, trace2]
layout = go.Layout(
    title="Number of Unknowns in the PPE System",
    yaxis=dict(
        title="Number of Unknowns"
    ),
    xaxis=dict(
        title="Frame"
    ),
    # legend=dict(orientation="h")
)

fig = go.Figure(data=data, layout=layout)
pio.write_image(fig, "ppe_unknowns.pdf")
py.plot(fig)