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
        return -1
    #returns lowest electronegativity atom
    central_atom = min(centre_candidates,key=lambda x: x[1])[0]
    return central_atom

#returns a list of substituents on central atom 
#from formula  "CH4" returns [H,H,H,H], "Fe(H20)6" returns [H20,H20,H20,H20,H20,H20]
def Substituents(formula):
    #create list of substituents ignoring central atom
    central_atom = CentralAtom(formula)
    remove_centre = formula.replace(central_atom,"",1)
    
    search = r"\((.*?)\)(\d*)|([A-Z][a-z]?)(\d*)"
    substituents_tuple = re.findall(search,remove_centre)

    substituents = []
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
