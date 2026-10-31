import matplotlib.pyplot as plt
import numpy as np
import VSEPR_function_v2 as vsepr

fig = plt.figure()
ax = fig.add_subplot(projection='3d')


#Sphere map
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

points,labels = vsepr.FindCoords(formula="HCN 2-")

#plot a sphere at each point generated if it has a bonded pair
for j,i in enumerate(points):
    if labels[j]!="Lone Pair":
        x = (10 * np.outer(np.cos(u), np.sin(v))) +i[0]*10
        y = (10 * np.outer(np.sin(u), np.sin(v))) +i[1]*10
        z = (10 * np.outer(np.ones(np.size(u)), np.cos(v))) +i[2]*10
        
        # Plot the surface
        ax.plot_surface(x, y, z)
        
#central atom sphere
x = (15 * np.outer(np.cos(u), np.sin(v))) 
y = (15 * np.outer(np.sin(u), np.sin(v))) 
z = (15 * np.outer(np.ones(np.size(u)), np.cos(v)))
        
# Plot the surface
ax.plot_surface(x, y, z)

# Set an equal aspect ratio
ax.set_aspect('equal')

plt.show()