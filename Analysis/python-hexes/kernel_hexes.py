import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import rgb2hex
import numpy as np
import pandas as pd
import math

# Makes hex grids for the dispersal kernel data
# Just change from calculating relative concentration to probability of finding all molecules in the same place at once

# Read data

data = pd.read_csv("Data/4-kernel/kernel_3_1_false_3739818.csv")

# Regenerate the coordinates

omega = end_rows.shape[0]-1 # -1 bc coordinate zero is global concentrations
R = int((1/6) * (math.sqrt((12*omega) - 3) - 3)) # rearranged equation to get the radius

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
if end_rows.shape[0]-1 != len(coord):
    print("Dimensions incorrect for this data file\n")
    print("Omega = "); print(omega)
    print("File has "); print(end_rows.shape[0]-1); print("pixels\n")



# Make the color map
gradient = ["#BEBEBE", "#005AB5"]
color_map = LinearSegmentedColormap.from_list("custom_gradient", gradient)
print(rgb2hex(color_map(0.5))) # Check

# Make colors based on occupation probabilities

colors = []
for x in range(len(coord)):
    #colors.append(["Gray"])
    current_row = end_rows[end_rows["pixel"] == x+1]
    da = "diff"

    if (current_row["sp_r1_1_"+da] >= 1) & (current_row["sp_r1_2_"+da] >= 1) & (current_row["sp_r1_3_"+da] >= 1):
        
    
    current_A = current_row["sp_r1_1_"+da] + current_row["sp_r1_2_"+da] + current_row["sp_r1_3_"+da]
    current_NA = current_row["sp_r2_1_"+da] + current_row["sp_r2_2_"+da] + current_row["sp_r2_3_"+da]


    current_A = current_A.to_numpy()[0]
    current_NA = current_NA.to_numpy()[0]

    current_rel_con = (current_A) / (current_A + current_NA)
    if np.isnan(current_rel_con):
        colors.append(["#BEBEBE"]) # empty space with nothink in it
    else:
        colors.append([rgb2hex(color_map(current_rel_con))]) # something from the color map

# Convert coordinates into something plottable
hcoord = [c[0] for c in coord]
vcoord = [2. * np.sin(np.radians(60)) * (c[1] - c[2]) /3. for c in coord]

# Make figure
fig, ax = plt.subplots(1)
ax.set_aspect('equal')

for x, y, c in zip(hcoord, vcoord, colors):
    color = c[0].lower()
    hex = RegularPolygon((x, y), numVertices=6, radius=2. / 3., 
                         orientation=np.radians(30), 
                         facecolor=color, alpha=1, edgecolor='k')
    ax.add_patch(hex)

# Make points so scope shows all the hexagons
ax.scatter(hcoord, vcoord, c=[c[0].lower() for c in colors], alpha=0)

plt.show()
                   