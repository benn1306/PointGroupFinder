from VSEPR_functions_v1 import groups,GroupNumFinder
import re
import numpy as np

#VSPER rules: 
#1. Find central atom - usually least electronegative
#2. Find number of bonding regions
#3. Find lone pairs
#4. position bonding regions and lone pairs as spaced out around central atom as possible
#5. calculate atom positions

#define repulsion constants for Thomson Problem
lp_repulsion = 1.25
bp_repulsion = 1.0

def ListElements(formula):
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

#Return most likely central atom
def CentralAtom(formula):
    elements = ListElements(formula)

    #For loop to search for potential central atoms
    central_cand = []
    for element in elements:
        if element[0] != "H" and not len(element) > 1: #Checks element is not H and not a substituent molecule
            central_cand.append([element,GroupNumFinder(element[0])])
    
    #Finds lowest electro negativity candidate and assume that is center
    central_atom = min(central_cand, key=lambda x: x[1])
    return central_atom

#Find number of bonding regions
def NumBondingRegions(formula):
    return len(ListElements(formula))-1

def NumLonePairs(formula):
    group = CentralAtom(formula)[1]
    used_ve = NumBondingRegions(formula)
    #Use group number to find valence electron count
    if 3<group<=12: return 0
    ve_count = group - 10 if group > 12 else group
    used_ve = NumBondingRegions(formula)
    return (ve_count + FindCharge(formula) - used_ve)//2

#Finds coordinates around sphere that minimise repulsion between charges (Thomson Problem)
def FindCoords(lp,bp,formula,steps=1000,lr=0.01):
    N = lp + bp
    if N < 2:
        raise ValueError("Not enough electron dense regions to find geometry")
    
    #array of repulsion values with length N and bp:lp of bp_repulsions and lp_repulsions
    charges = np.array([bp_repulsion] * bp + [lp_repulsion] * lp)
    clean_formula = re.sub(CentralAtom(formula)[0][0],"",formula)
    labels = ["".join(i)+" "+str(x+1) for x,i in enumerate(ListElements(clean_formula))] + ["Lone Pair"] * lp  #format numbered labels for bonded atoms

    #Start at random positions
    coords = np.random.normal(size=(N,3))
    coords /= np.linalg.norm(coords, axis=1, keepdims=True) #Force points onto surface of sphere

    #Iterative relaxation of repulsion
    for _ in range(steps):
        forces = np.zeros_like(coords)
        for i in range(N):
            for j in range(N):
                if i == j:
                    continue #dont check repulsion of the same coordinate

                #vector pushes i away from j
                diff = coords[i] - coords[j]
                dist = np.linalg.norm(diff)

                if dist > 1e-5:
                    #Use variant of coulombs law to find force on coordinate i
                    forces[i] += (charges[i]*charges[j]*diff/(dist**3))

        #move along net repulsive force vector            
        coords += lr * forces
        coords /= np.linalg.norm(coords,axis=1,keepdims=True) #ensure everything is still on sphere surface

    return coords, labels


#prints list of angles between bonds (work on)
def print_bond_angles(formula):
    coords,labels = FindCoords(NumLonePairs(formula),NumBondingRegions(formula),formula) #generate spaced out coordinates and associated labels for given formula
    N = len(coords) #How many bond angles need to be printed

    #Formating
    print(f"\n{'Pair I':<15} {'Pair J':<15} {'Angle (Degrees)':<15}")
    print("-" * 48)

    #Iterate through coords comparing dot products to work out angle in degrees and print in format
    for i in range(N):
        for j in range(i + 1, N):
            dot_product = np.dot(coords[i], coords[j])
            # Clip to avoid float precision issues with acccos
            angle = np.degrees(np.arccos(np.clip(dot_product, -1.0, 1.0)))
            print(f"{labels[i]:<15} {labels[j]:<15} {angle:<15.2f}")



if __name__ == "__main__":
    print_bond_angles(input("Enter formula to generate bond angles"))