import re
#Dictionary of most elements and their group number
groups = {1:["H","Li","Na","K","Rb","Cs","Fr"]
          ,2:["Be","Mg","Ca","Sr","Ba","Ra"]
          ,3:["Sc","Y"],4:["Ti","Zr","Hf"]
          ,5:["V","Nb","Ta"],6:["Cr","Mo","W"]
          ,7:["Mn","Tc","Re"],8:["Fe","Ru","Os"]
          ,9:["Co","Rh","Ir"],10:["Ni","Pd","Pt"]
          ,11:["Cu","Ag","Au"],12:["Zn","Cd","Hg"]
          ,13:["B","Al","Ga"]
          ,14:["C","Si","Ge","Sn"]
          ,15:["N","P","As","Sb","Bi"]
          ,16:["O","S","Se","Te","Po","Lv"]
          ,17:["F","Cl","Br","I","At","Uus"]
          ,18:["He","Ne","Ar","Kr","Xe","Rn","Uuo"]}

#Dictionary of shapes, keys given (lone pairs, bonding pairs)
shapes = {(0, 2): "Linear",(3, 2): "Linear",(4, 2): "Linear",
    (1, 2): "Bent",(2, 2): "Bent",
    (0, 3): "Trigonal Planar",(1, 3): "Trigonal Pyramid",
    (2, 3): "T-Shape",(3, 3): "T-Shape",
    (0, 4): "Tetrahedral",
    (1, 4): "Sawhorse",
    (2, 4): "Square Planar",
    (0, 5): "Trigonal Bipyramidal",
    (1, 5): "Square Pyramid",
    (0, 6): "Octahedral"}

#Find out if a substituent will form single, double or triple bonds
def BondOrderApprox(substituent):
    group_num = GroupNumFinder(substituent)
    if group_num <=12 or len(substituent)>2:
        return 1
    else:
        num_ve = group_num - 10
        return 8-num_ve

#Finds group number of given element (if not in dict return -1)
def GroupNumFinder(element):
    for group_num, elements in groups.items():
        if element in elements:
            return group_num
    return -1

#Find central atom
def CentralAtom(formula):
    #create list of elements, ignoring parts in brackets
    remove_ligands = re.sub(r"\(.*?\)\d*", "", formula)
    elements = re.findall(r"[A-Z][a-z]?",remove_ligands)
    
    #Define list to put 
    centre_candidates = []

    for element in elements:
        #Central atom can't be H and group 18 functionality to be added later
        if element == "H":
            continue
        group_num = GroupNumFinder(element)
        if group_num == 18:
            continue

        centre_candidates.append((element,group_num))
    #if no candidates return -1 to indicate error
    if not centre_candidates:
        return ""
    #returns lowest electronegativity atom
    central_atom = min(centre_candidates,key=lambda x: x[1])[0]
    return central_atom

#returns a list of substituents on central atom 
#from formula  "CH4" returns [H,H,H,H], "Fe(H20)6" returns [H20,H20,H20,H20,H20,H20]
def Substituents(formula):
    substituents = []
    #create string of substituents ignoring central atom
    central_atom = CentralAtom(formula)
    central_atom_count = re.search(rf"{central_atom}(\d*)",formula)
    try:
        if central_atom_count:
            substituents.extend([central_atom]*(int(central_atom_count.group(1))-1))
    except:
        pass
    remove_centre = formula.replace(central_atom,"",1)
    
    #Finds all substituents in the compound and how many times they occur
    search = r"\((.*?)\)(\d*)|([A-Z][a-z]?)(\d*)"
    substituents_tuple = re.findall(search,remove_centre)

    #Unpacks encoded substituents
    for molecule,molecule_quant,atom,atom_quant in substituents_tuple:
        if molecule:
            if molecule_quant:
                num = int(molecule_quant)  
            else:num = 1
            
            substituents.extend([molecule]*num)

        elif atom:
            if atom_quant:
                num = int(atom_quant)
            else:num = 1

            substituents.extend([atom]*num)

    return substituents

#Finds charge of formula
def ChargeFinder(formula):
    charge_search = re.search(r"\]?(\d*)([+-])$",formula.strip())
    if not charge_search:
        return 0
    
    charge, sign = charge_search.groups()
    if charge:
        charge = (int(charge))
    else:
        charge = 1
    if sign == "-":
        charge *= -1
    return charge

#Finds number of lone pairs
def LonePairsBondedPairs(formula):
    
    subs = Substituents(formula)
    num_subs = 0
    for sub in subs:
        num_subs += BondOrderApprox(sub)

    #group number
    central_atom = CentralAtom(formula)
    num_ve = GroupNumFinder(central_atom)

    #returns molecule charge
    charge_e = ChargeFinder(formula)

    #Exceptions for transition metals
    if 3<= num_ve <= 12:
        if num_subs == 4:
            if central_atom in groups[num_ve][1:]:
                return(2,4) 
            #Will add ligand evaluation to allow for period 4 transition metals to have accurate shapes
        return(0,num_subs)

    #Valence electrons
    if num_ve > 12:
        num_ve -= 10

    num_ve -= charge_e
    #Number of lone pairs based on charge + valence electrons - substituents
    num_lp = (num_ve - num_subs)//2
    return (max(-1,num_lp),num_subs)

def FindShape(formula):
    lp_bp = LonePairsBondedPairs(formula)
    if lp_bp[0] == -1:
        return "Formula Error"
    print(lp_bp)
    return shapes.get(lp_bp,"Shape Not Defined")

if __name__ == "__main__":
    run = "y"
    while run != "n":
        print(FindShape(str(input("Enter simple molecule: "))))
        run = input("Continue (y/n): ")
