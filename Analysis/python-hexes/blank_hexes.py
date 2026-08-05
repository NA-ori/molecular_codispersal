import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import rgb2hex
import numpy as np
import pandas as pd
import math

# Blank hexes, for fun

# Read data

# Regenerate the coordinates

D = 3
s = 2
R = s*D + (D//2)
omega =  omega = 3*R*(R+1) + 1
R = int((1/6) * (math.sqrt((12*omega) - 3) - 3)) # rearranged equation to get the radius

# Same as matlab coord generation
coord = []
for x in range(-R, R+1):
    for y in range(-R, R+1):
        for z in range(-R, R+1):
            if (x+y+z) == 0:
                coord.append([x,y,z])


# Make colors based on relative concentrations

da = "ad"     # diff or ad, whichever you want to graph

colors = []
for x in range(len(coord)):
    colors.append(["#BEBEBE"]) # empty space with nothink in it


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
                   