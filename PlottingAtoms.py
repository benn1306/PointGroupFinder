import matplotlib.pyplot as plt
import numpy as np
import VSEPR_function_v2 as vsepr
import re

def GetVariables(element):
    """Returns a colour and size depending on group number of element"""
    group = vsepr.GroupNumFinder(element)
    if group == 1:
        colour = "white"
    elif group < 3:
        colour = "purple"
    elif group <= 13:
        colour = "#DAA48F" #tan
    elif group == 14:
        colour = "#1B141A" #charcoal
    elif group == 15:
        colour = "blue"
    elif group == 16:
        colour = "red"
    elif group == 17:
        colour = "green"
    else:
        colour = "yellow"
    
    size = 0.4
    return size, colour

def PlotAtoms(formula):
    """Plots spheres for each point an atom occurs"""
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    #Sphere map
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)

    #Generate all coordinates scaled around central atom at (0,0,0)
    points, labels = vsepr.FindCoords(formula)

    #Dictionary to track position of subgroup center
    subgroup_centers = {}

    #Locates subgroup centers coordinates
    for j, label in enumerate(labels):
        if " -> " in label:
            parts = label.split(" -> ")
            group_name = parts[0] #Full subgroup
            atom_type = parts[1] 
            if not re.search(r"\s+\d+", atom_type) and atom_type != "Lone Pair":
                subgroup_centers[group_name] = points[j]

    #Plots all atoms and connects them with lines
    for j, i in enumerate(points):
        if labels[j] == "Lone Pair": #Dont draw lone pairs
            continue
            
        cx, cy, cz = i[0], i[1], i[2]
        
        #Extracting element from label
        if " -> " in labels[j]:
            element = labels[j].split(" -> ")[1]
            if re.search(r"\d+", element):
                element = "".join(element.split()[:-1])
        else:
            element = labels[j].split(" ")[0]
        
        #Draw atom sphere at appropriate location
        r,colour = GetVariables(element)
        x = (r * np.outer(np.cos(u), np.sin(v))) + cx
        y = (r * np.outer(np.sin(u), np.sin(v))) + cy
        z = (r * np.outer(np.ones(np.size(u)), np.cos(v))) + cz
        ax.plot_surface(x, y, z, alpha=0.8, color=colour)
        
        #Label atom with appropriate label
        ax.text(cx * 1.1, cy * 1.1, cz * 1.1, labels[j], color="black", fontsize=8)

        #Lines between atoms (bonds)
        if " -> " in labels[j]:
            parts = labels[j].split(" -> ")
            group_name = parts[0]
            atom_type = parts[1]
            
            #If looking at a subgroup center, connect to central atom and dont bother connecting to self
            if not re.search(r"\s+\d+", atom_type) and atom_type != "Lone Pair":
                ax.plot3D([0, cx], [0, cy], [0, cz], color="gray", linewidth=0.5)
            else:
                #Connect atoms in subgroups to center of subgroup
                if group_name in subgroup_centers:
                    parent_xyz = subgroup_centers[group_name]
                    ax.plot3D([parent_xyz[0], cx], [parent_xyz[1], cy], [parent_xyz[2], cz], color="black", linewidth=0.5)
        else:
            #Connect standalone atoms to central atom
            ax.plot3D([0, cx], [0, cy], [0, cz], color="gray", linewidth=0.5)

    #Plot the central atom at (0,0,0)
    central_atom = vsepr.CentralAtom(formula)[0][0]
    r,colour = GetVariables(central_atom)
    x_c = r * np.outer(np.cos(u), np.sin(v))
    y_c = r * np.outer(np.sin(u), np.sin(v))
    z_c = r * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x_c, y_c, z_c, color=colour)
    ax.text(0, 0, 0.9, central_atom, color="black", weight="bold", ha="center")

    #Label axis and set aspect
    ax.set_aspect("equal")
    ax.set_xlabel("X Axis")
    ax.set_ylabel("Y Axis")
    ax.set_zlabel("Z Axis")

    plt.show()


if __name__ == "__main__":
    PlotAtoms(input("enter molecule to display: "))
