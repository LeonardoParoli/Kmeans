import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the file and extract the coordinates
with open("initial_points.txt", "r") as file:
    lines = file.readlines()
    points = [eval(line.strip()) for line in lines]

# Separate the x, y, z coordinates
x_coords = [point[0] for point in points]
y_coords = [point[1] for point in points]
z_coords = [point[2] for point in points]

with open("real_centroids.txt", "r") as file:
    lines = file.readlines()
    additional_points = [eval(line.strip()) for line in lines]

# Separate the x, y, z coordinates of additional points
additional_x_coords = [point[0] for point in additional_points]
additional_y_coords = [point[1] for point in additional_points]
additional_z_coords = [point[2] for point in additional_points]
# Plot the points in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_coords, y_coords, z_coords, color='blue', s=0.05)
ax.scatter(additional_x_coords, additional_y_coords, additional_z_coords, color='red', s=25)


# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Scatter Plot')

# Show the plot
plt.show()