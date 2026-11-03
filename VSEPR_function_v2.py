#VSPER rules: 
#1. Find central atom - usually least electronegative
#2. Find number of bonding regions
#3. Find lone pairs
#4. position bonding regions and lone pairs as spaced out around central atom as possible
#5. calculate atom positions
import re
import numpy as np

#Dictionary of most elements and their group number
groups = {1:["H","Li","Na","K","Rb","Cs","Fr"]
          ,2:["","Be","Mg","Ca","Sr","Ba","Ra"]
          ,3:["","","","Sc","Y"],4:["","","","Ti","Zr","Hf"]
          ,5:["","","","V","Nb","Ta"],6:["","","","Cr","Mo","W"]
          ,7:["","","","Mn","Tc","Re"],8:["","","","Fe","Ru","Os"]
          ,9:["","","","Co","Rh","Ir"],10:["","","","Ni","Pd","Pt"]
          ,11:["","","","Cu","Ag","Au"],12:["","","","Zn","Cd","Hg"]
          ,13:["","B","Al","Ga"]
          ,14:["","C","Si","Ge","Sn"]
          ,15:["","N","P","As","Sb","Bi"]
          ,16:["","O","S","Se","Te","Po","Lv"]
          ,17:["","F","Cl","Br","I","At","Uus"]
          ,18:["He","Ne","Ar","Kr","Xe","Rn","Uuo"]}

#define repulsion constants for Thomson Problem
lp_repulsion = 1.25
bp_repulsion = 1.0

def GroupNumFinder(element):
    """Returns the group number of given element, -1 if not found"""
    for group_num, elements in groups.items():
        if element in elements:
            return group_num
    return -1


def ListElements(formula):
    """Breaks a formula into sublists of elements for each substituent"""
    #list of all elements and substituents encoded
    element_encoded = re.findall(r"[A-Z][a-z]?\d*|\(.*?\)\d*", formula)
    
    #for loop to expand terms
    elements = []
    for sub in element_encoded:
        #Expanding bracketed groups
        if sub.startswith("("):
            match = re.match(r"\((.*?)\)(\d*)", sub)
            sub_formula = match.group(1) #Inside brackets

            #int outside brackets
            if match.group(2):
                sub_quant = int(match.group(2))  
            else:
                sub_quant = 1  
            
            #Expanding substituent within intself
            inner_atoms = re.findall(r"[A-Z][a-z]?\d*", sub_formula)
            sub_list = []
            for atom in inner_atoms:
                #Find how many of each atom
                count_atom = re.search(r"\d+", atom) 
                if count_atom:
                    count = int(count_atom.group())
                else:
                    count = 1

                #Find element alone and add to the list count times
                element = re.sub(r"\d+", "", atom)
                sub_list.extend([element] * count)
            
            #Append sub list to elements sub_quant times
            for _ in range(sub_quant):
                elements.append(sub_list)
            
        #Expanding stand alone atoms
        else:
            #How many of each independent atom
            count_match = re.search(r"\d+", sub)
            if count_match:
                count = int(count_match.group()) 
            else:
                count = 1

            #Element without number
            element = re.sub(r"\d+", "", sub)
            
            #Put element in independent list and add to overall list count times
            for _ in range(count):
                elements.append([element])
            
    return elements

def FindCharge(formula):
    """Returns the charge of a molecule"""
    #check if final number followed by +/- 
    charge_search = re.search(r"\]?(\d*)([+-])$",formula.strip())
    if not charge_search: #If no charge given return 0
        return 0
    
    #assign the magnitude and sign of the charge
    charge, sign = charge_search.groups()
    if charge:
        charge = (int(charge))
    else:
        charge = 1
    #Multiply by correct sign and return as int
    if sign == "-":
        charge *= -1
    return charge

def CentralAtom(formula):
    """Returns most likely central atom of a molecule based on electronegativity"""
    elements = ListElements(formula)

    #For loop to search for potential central atoms
    central_cand = []
    for element in elements:
        if element[0] != "H" and not len(element) > 1: #Checks element is not H and not a substituent molecule
            central_cand.append([element,GroupNumFinder(element[0])])
    
    #Finds lowest electro negativity candidate and assume that is center
    central_atom = min(central_cand, key=lambda x: x[1])
    return central_atom

def NumBondingRegions(formula):
    """Returns number of bonded regions about a central atom in a molecule"""
    return len(ListElements(formula))-1

def NumLonePairs(formula):
    """Returns the number of lone pairs on the central atom of a molecule"""
    group = CentralAtom(formula)[1]
    used_ve = NumBondingRegions(formula)

    #Use group number to find valence electron count
    if 3<group<=12: return 0
    ve_count = group - 10 if group > 12 else group

    remaining_e = ve_count + FindCharge(formula) - used_ve

    #round up radicals to be equivalent to lone pair
    return remaining_e // 2 if remaining_e%2==0 else (remaining_e //2) + 1


def GetRotationMatrix(vec_source, vec_target):
    """Returns rotation matrix for directing subgroups away from central atom"""
    #calculate length of target vector and swing direction
    norm_s = np.linalg.norm(vec_source)
    norm_t = np.linalg.norm(vec_target)
    if norm_s < 1e-5 or norm_t < 1e-5: #if either length is less than 1e-5 return identity matrix
        return np.eye(3)
        
    #Normalize both vectors to only look at rotation angle
    u = vec_source / norm_s
    v_target = vec_target / norm_t
    
    #calculate dot product (cos of angle between the vectors)
    dot = np.dot(u, v_target)
    if dot > 0.9999: return np.eye(3) #if pointing same direction not rotation
    if dot < -0.9999: return -np.eye(3) #if opposite direction 180 flip
    
    v = np.cross(u, v_target) #axis of rotation
    s = np.linalg.norm(v) #length of v is sin of angle between them
    
    #Skew-symmetric matrix with components of v
    kmat = np.array([
        [  0,   -v[2],  v[1]],
        [ v[2],   0,   -v[0]],
        [-v[1],  v[0],   0  ]
    ])
    
    #Rodigues rotation formula R = I + K + K^2(1-cos / sin^2)
    return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - dot) / (s ** 2))


def HandleSubGroups(coords, labels, formula):
    """Returns coordinates for elements in subgroups"""
    #empty lists to update coordinates without using excess memory by appending to np array
    final_coords = []
    final_labels = []
    
    #iterate through list of 
    for i, label in enumerate(labels):
        if len(label) > 3 and not label == "Lone Pair" and " -> " not in label:
            clean_sub_formula = "".join(label.split())[:-1] #removes number from end of substituent
    
            #Stores the label for the group including the number of times it had previously occured
            unique_group_id = f"{label}"
            
            #Find coordinates of substituent elements about point (0,0,0)
            new_coords, new_labels = FindCoords(clean_sub_formula)
            
            #Defines the location of the central atom of the substituent
            parent_bond_vector = coords[i]
            global_ligand_center = parent_bond_vector * 2.0  
            
            #Find the center axis of whatever outer atoms exist and point them outward
            sub_bp = NumBondingRegions(clean_sub_formula) + 1 if not 3 < CentralAtom(formula)[1] <= 12 else NumBondingRegions(clean_sub_formula) 
            local_axis_source = np.mean(new_coords[:sub_bp], axis=0)
            target_axis_heading = parent_bond_vector
            
            #Generate matrix to calculate rotation
            rot_matrix = GetRotationMatrix(local_axis_source, target_axis_heading)
            
            #Finds central atom of substituent
            sub_center_info = CentralAtom(clean_sub_formula)
            sub_center_atom = sub_center_info[0][0] 
            
            #append coordinate and label for centre of subgroup
            final_coords.append(global_ligand_center)
            final_labels.append(f"{unique_group_id} -> {sub_center_atom}")
            
            #append new coordinates for subgroup elements rotated to point away from main central atom
            for (s_coord, s_label) in zip(new_coords, new_labels):
                if s_label == "Lone Pair": 
                    continue
                global_xyz = rot_matrix.dot(s_coord * 0.96) + global_ligand_center
                final_coords.append(global_xyz)
                final_labels.append(f"{unique_group_id} -> {s_label}")
        else: #if not a subgroup just append given label and coordinate
            final_coords.append(coords[i] * 2.0)
            final_labels.append(label)
    
    #convert to np array and return
    return np.array(final_coords), final_labels


def FindCoords(formula, steps=1000, lr=0.01):
    """Finds coordinates around sphere that minimise repulsion between charges (Thomson Problem)"""
    lp, bp = NumLonePairs(formula), NumBondingRegions(formula)
    N = lp + bp
    
    #Array of length N consisting of the relevant like charges for the approximation
    charges = np.array([bp_repulsion] * bp + [lp_repulsion] * lp)
    
    centeral_atom = CentralAtom(formula)[0][0]
    clean_formula = re.sub(re.escape(centeral_atom), "", formula, count=1) #Remove first instance of central atom from formula
    
    #List names of substituents and lone pairs
    labels = ["".join(i)+" "+str(x+1) for x,i in enumerate(ListElements(clean_formula))] + ["Lone Pair"] * lp  

    #Generate N random starting point coordinates on the surface of a sphere
    coords = np.random.normal(size=(N,3))
    coords /= np.linalg.norm(coords, axis=1, keepdims=True) 

    #Iterative relaxation for set amount of steps
    for _ in range(steps):
        forces = np.zeros_like(coords) #array of zeros that has same structure as coords array
        for i in range(N):
            for j in range(N):
                if i == j: continue #Dont compare coordinates with themselves
                #Find differentce in coordinates and convert to distance
                diff = coords[i] - coords[j] 
                dist = np.linalg.norm(diff)
                #Vector that they i will be repelled from j by
                if dist > 1e-5: 
                    forces[i] += (charges[i]*charges[j]*diff/(dist**3))
         
        coords += lr * forces #adds the force vectors to the random coordinates
        coords /= np.linalg.norm(coords,axis=1,keepdims=True) #Snaps coordinates to sphere surface

    coords, labels = HandleSubGroups(coords, labels, formula)
    return coords, labels

if __name__ == "__main__":
    print(FindCoords(input("Enter formula to generate bond angles: ")))