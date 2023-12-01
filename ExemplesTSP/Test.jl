using Plots

# Coordinates of the central point
center_point = (0, 0)
xlist = [1,2,3]
ylist = [1,2,3]

# Coordinates of other points
other_points = [(1, 1), (2, 2), (3, 1), (1, -2)]

# Extract x and y coordinates for plotting
center_x, center_y = center_point
other_x, other_y = zip(other_points...) |> collect
#scatter([center_x], [center_y], color=:red, label="Center Point")
println(other_x)
scatter(other_x, other_y, color=:blue, label="Other Points")

# Plot the central point
# scatter!(p,[center_x], [center_y], color=:red, label="Center Point")

# Plot the other points
#scatter!(p,other_x, other_y, color=:blue, label="Other Points")

# Connect the central point to other points with lines
# for point in other_points
#     plot!(p,[center_x, point[1]], [center_y, point[2]], line=:solid, color=:gray)
# end

# Add labels and legend
# title!("Connecting Points")
# xlabel!("X-axis")
# ylabel!("Y-axis")

# Show the plot
# display(p)