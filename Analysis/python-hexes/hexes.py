import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import rgb2hex
import numpy as np
import pandas as pd

# Some settings based on the world size
# This will make sure coordinates are accurate
# Manually input this before running

separation_distance = 7
sites = 2

R = separation_distance*sites + (separation_distance // 2) # Calculate world radius
omega = 3*R*(R+1) + 1

# Same as matlab coord generation
coord = []
for x in range(-R, R+1):
    for y in range(-R, R+1):
        for z in range(-R, R+1):
            if (x+y+z) == 0:
                coord.append([x,y,z])

# Checks
print(R)
print(omega)
print(len(coord)) # should be the same as omega

# Read data

data = pd.read_csv("Data/raw/2-comp/curves_s2/complete_curve_4829.csv")
end_rows = data[data['t'] > 3499] # Just in case there's some variation in the exact final time

# Checks
if end_rows.shape[0]-1 != len(coord):
    print("Dimensions incorrect for this data file\n")
    print("Omega = "); print(omega)
    print("File has "); print(end_rows.shape[0]-1); print("pixels\n")


# Make the color map
gradient = ["#D62800", "#1AB3FF"]
color_map = LinearSegmentedColormap.from_list("custom_gradient", gradient)
print(rgb2hex(color_map(0.5))) # Check

# Make colors based on relative concentrations
colors = []
for x in range(len(coord)):
    #colors.append(["Gray"])
    current_row = end_rows[end_rows["pixel"] == x+1]
    # Look at adsorbed
    current_A = current_row["sp_r1_1_ad"] + current_row["sp_r1_2_ad"] + current_row["sp_r1_3_ad"]
    current_NA = current_row["sp_r2_1_ad"] + current_row["sp_r2_2_ad"] + current_row["sp_r2_3_ad"]
    # Look at diffused
    current_A = current_row["sp_r1_1_diff"] + current_row["sp_r1_2_diff"] + current_row["sp_r1_3_diff"]
    current_NA = current_row["sp_r2_1_diff"] + current_row["sp_r2_2_diff"] + current_row["sp_r2_3_diff"]

    current_A = current_A.to_numpy()[0]
    current_NA = current_NA.to_numpy()[0]

    try:
        current_rel_con = (current_A) / (current_A + current_NA)
    except ZeroDivisionError:
        current_rel_con = 512 # something nonsensical, also a reference :3
    if current_rel_con != 512:
        colors.append([rgb2hex(color_map(current_rel_con))])
    elif current_rel_con == 512:
        colors.append(["Gray"])


# Convert coordinates into something plottable
hcoord = [c[0] for c in coord]
vcoord = [2. * np.sin(np.radians(60)) * (c[1] - c[2]) /3. for c in coord]

# Produce figure
fig, ax = plt.subplots(1)
ax.set_aspect('equal')

for x, y, c in zip(hcoord, vcoord, colors):
    color = c[0].lower()
    hex = RegularPolygon((x, y), numVertices=6, radius=2. / 3., 
                         orientation=np.radians(30), 
                         facecolor=color, alpha=0.2, edgecolor='k')
    ax.add_patch(hex)

# Also add scatter points in hexagon centres
ax.scatter(hcoord, vcoord, c=[c[0].lower() for c in colors], alpha=0)

plt.show()
                   