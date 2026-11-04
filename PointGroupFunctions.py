import numpy as np
import VSEPR_function_v2 as vsepr

class Molecule:
    def __init__(self,formula):
        """initialise object to find all symmetry functions"""
        self.formula = formula
        self.coords,self.labels = vsepr.FindCoords(self.formula,sub_group_coords=False)

        #remove lone pair coords from labels and coords
        del_lone_pair_idx = [i for i,label in enumerate(self.labels) if label != "Lone Pair"]
        self.coords = [self.coords[i] for i in del_lone_pair_idx]
        self.labels = ["".join(self.labels[i].split()[:-1]) for i in del_lone_pair_idx]

    def MirrorPlanes(self,tol=1e-4):
        """Returns list of mirror planes as vectors"""
        candidates = []
        if not self.coords:
            return []
        
        for i in range(0,len(self.coords)):
            for j in range(i+1,len(self.coords)):
                #find normal vectors that cut through 2 points and the origin
                vector = np.cross(self.coords[i],self.coords[j])
                normalised_vector = np.linalg.norm(vector)
                if normalised_vector > 0.1: #0.1 not tol to filter out noise in opposite atom cross vectors
                    candidates.append(vector / normalised_vector)

                #find normal vectors that cut through two identical atoms 
                if self.labels[i] == self.labels[j]:
                    vector = self.coords[i]-self.coords[j]
                    normalised_vector = np.linalg.norm(vector)
                    if normalised_vector > tol:
                        candidates.append(vector / normalised_vector)

        unique_candidates = []
        #remove repeat values   
        for i in candidates:
            i_standard = np.copy(i)

            #Find the index of the component with the largest absolute value
            max_idx = np.argmax(np.abs(i_standard))
            
            #Standardise based on its sign
            if i_standard[max_idx] < 0:
                i_standard = -i_standard
                
            #slightly relaxed tolerance to prevent vsepr coordinate drift allowing same mirror planes through twice
            if not any(np.allclose(i_standard, u, atol=0.05) for u in unique_candidates):
                unique_candidates.append(i_standard)

        #validating reflections
        validated_candidates = []
        for cand in unique_candidates:
            R = np.eye(3) - 2 * np.outer(cand,cand)
            coords_arr = np.array(self.coords)
            reflected_coords = coords_arr.dot(R.T)
            is_mirror = True
            for ref_coord, ref_label in zip(reflected_coords, self.labels):
                #Calculate distances from this single reflected atom to all original atoms
                distances = np.linalg.norm(coords_arr - ref_coord, axis=1)
                
                #Find indices where the distance is small and the label matches
                matches = [idx for idx, dist in enumerate(distances) if dist < 0.05 and self.labels[idx] == ref_label]
                
                #If any of the new coordinates have no matching partner in the original coords, mirror plane is invalid
                if not matches:
                    is_mirror = False
                    break
                    
            if is_mirror:
                validated_candidates.append(cand)
        
        return validated_candidates
    
    def RotationalAxis(self):
        """returns list of rotational axis"""
        return
    def PrincipalAxis(self):
        """returns highest rotational axis"""
        return
    def Dihedral(self):
        """returns dihedral C2 axis"""
        return
    def InversionCenter(self):
        """returns True/False depending on if there is an inversion center"""
        return
    def HighSymmetry(self):
        """Returns True and symmetry label if high symmetry else returns False"""
        return
    
mol = Molecule("BH4")
print(len(mol.MirrorPlanes()))