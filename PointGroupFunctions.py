import numpy as np
import VSEPR_function_v2 as vsepr

class Molecule:
    def __init__(self,formula):
        """initialise object to find all symmetry functions"""
        #define formula and all coords and labels including lone pairs
        self.formula = formula
        self.coords,self.labels = vsepr.FindCoords(self.formula,sub_group_coords=False)

        #remove lone pair coords from labels and coords
        del_lone_pair_idx = [i for i,label in enumerate(self.labels) if label != "Lone Pair"]
        self.coords = [self.coords[i] for i in del_lone_pair_idx]
        self.coords_arr = np.array(self.coords)
        self.labels = ["".join(self.labels[i].split()[:-1]) for i in del_lone_pair_idx]
        self.central_atom = vsepr.CentralAtom(self.formula)[0][0]

        #define valid mirror planes and rotational axes and principal axis
        self.mirrors = self.MirrorPlanes()
        self.rotational_axes = self.RotationAxes()
        self.principal = self.PrincipalAxis()

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
            reflected_coords = self.coords_arr.dot(R.T)            
            is_mirror = self.ComparePoints(reflected_coords)       
            if is_mirror:
                validated_candidates.append(cand)
        
        return validated_candidates
    
    def PrincipalAxis(self):
        """returns highest rotational axis"""
        if self.rotational_axes:
            return max(self.rotational_axes,key=lambda x:x[1])
        return 
    
    def DihedralAxes(self, tol=0.05):
        """returns list of dihedral C2 axes"""
        dihedral_axes = []
        for axis in self.rotational_axes:
            if not axis[1] == 2: #if not a C2 axis check next
                pass
            
            dot_prod = np.dot(self.principal[0],axis[0]) #if perpendicular this = ~0
            if dot_prod < tol:
                dihedral_axes.append(axis)

        return dihedral_axes 
    
    def InversionCenter(self):
        """returns True/False depending on if there is an inversion center"""
        if not self.coords:
            return False
        #flip all points and then compare new points with original array
        inverted_coords = -self.coords_arr
        return self.ComparePoints(inverted_coords)
    
    def HighSymmetry(self):
        """Returns True and symmetry label if high symmetry else returns False"""
        if len(self.labels) == 1: #check if only 2 atoms if so it is linear so high symmetry
            if self.labels[0] != self.central_atom:
                return True, "C inf,v"
            return True, "D inf,h"
        
        parts = set(self.labels)
        length = len(self.labels)

        if len(parts) == 1: #check if there is only one type of substituent, if so see if lone pairs
            if vsepr.NumLonePairs(self.formula) == 0:
                if length == 4: #no lone pairs and 4 of same atom = tetrahedral
                    return True, "T d"
                if length == 6: #no lone pairs and 6 of same atom = octahedral
                    return True, "O h"
        return False, "Not High Sym."
    
    def ComparePoints(self,new_coords):
        """Compare distances between points in a set of coords and original coords"""
        for ref_coord, ref_label in zip(new_coords, self.labels):
                #Calculate distances from this single reflected atom to all original atoms
                distances = np.linalg.norm(self.coords_arr - ref_coord, axis=1)
                
                #Find indices where the distance is small and the label matches
                matches = [idx for idx, dist in enumerate(distances) if dist < 0.05 and self.labels[idx] == ref_label]
                
                #If any of the new coordinates have no matching partner in the original coords, mirror plane is invalid
                if not matches:
                    return False
        return True
    
    def IntersectionAxes(self,tol=1e-2):
        """returns list of intersection axes of mirror planes"""
        if len(self.mirrors) < 2:
            return []
        
        axes = []
        for i in range(len(self.mirrors)):
            for j in range(i+1, len(self.mirrors)):
                #take cross product to check for intersection if zero planes are parallel
                axis = np.cross(self.mirrors[i],self.mirrors[j])
                axis_norm = np.linalg.norm(axis)
                if axis_norm > tol:
                    unit_axis = axis / axis_norm
                    
                    #Find the index of the component with the largest absolute value
                    max_idx = np.argmax(np.abs(unit_axis))
                    #Standardise based on its sign
                    if unit_axis[max_idx] < 0:
                        unit_axis = -unit_axis

                    #only add if its a new axis
                    if not any(np.allclose(unit_axis, u, atol=tol) for u in axes):
                        axes.append(unit_axis)
        return axes
    def IsHorizontalMirror(self,tol=0.05):
        """checks if mirror is perpendicular to principal rotation axis and returns True/False"""
        for mirror in self.mirrors:
            dot_prod = np.dot(mirror,self.principal[0])
            if dot_prod > tol:
                return True
        return False

    def RotationMatrix(self, axis, angle_degrees):
        """Returns a rotation matrix for a rotation around a given axis vector by an angle in degrees."""
        #Convert the angle to radians and calculate sin / cos
        theta = np.radians(angle_degrees)
        cos_t = np.cos(theta)
        sin_t = np.sin(theta)
        
        #skew-symmetric matrix with componenets of the axis
        kmat = np.array([
            [0,  -axis[2],   axis[1]],
            [ axis[2],     0,  -axis[0]],
            [-axis[1],   axis[0],     0]
        ])
        
        #Rodigues rotation formula R = I + K + K^2(1-cos / sin^2)
        return np.eye(3) + sin_t * kmat + (1 - cos_t) * np.dot(kmat, kmat)
    
    def RotationAxes(self, tol=0.05):
        """returns rotation axes and order of angle"""
        cand_axes = self.IntersectionAxes()
        verified_axes = []
        #check each axis for each order of rotation until it finds one or reaches gets past 2
        for axis in cand_axes:
            for order in [6,5,4,3,2]:
                angle = 360 // order

                #generate rotation matrix and rotated coordinates
                R = self.RotationMatrix(axis,angle)
                rotated_coords = self.coords_arr.dot(R.T)

                #append if rotated coords matches original coords
                if self.ComparePoints(rotated_coords):
                    verified_axes.append((axis,order))
                    break

        return verified_axes
            
    def FindPointGroup(self):
        """Goes through flow chart steps to determine point group"""
        if self.HighSymmetry()[0]: #check for high symmetry molecule
            return self.HighSymmetry()[1]
        
        if not self.rotational_axes: #if no rotational axes (E) check for inversion center
            if self.InversionCenter():
                return "C i" #Point group for symmetry E with inversion center
            return "C 1" #point group for symmetry E no inversion center no improper rotation axes
        
        #Find mirror label based on if there is a horizontal mirror plane
        mirror_label = ""
        if self.mirrors:
            mirror_label = "h" if self.IsHorizontalMirror() else "v"
        
        #find rotation label based on if there is a dihedral rotational axis
        rot_label = "D" if self.DihedralAxes() else "C"

        #put it all together and return
        sym_components = [rot_label,str(self.principal[1]),mirror_label]
        return " ".join(sym_components)
        
                
#test values and outputs to ensure correct point group being found
if __name__ == "__main__":
    test_molecules = ["SF6", "SF4Cl2", "PF5", "BF3", "SF5Cl", "CF3Cl", "CFClBrI","IF7","SF4"]
    
    print(f"{'Formula':<10} | {'Planes':<6} | {'Inv':<5} | {'Principal':<10} | {'Group Label'}")
    print("-" * 70)
    for formula in test_molecules:
        try:
            mol = Molecule(formula)
            planes = len(mol.MirrorPlanes())
            inv = mol.InversionCenter()
            
            p_axis = mol.PrincipalAxis()
            p_str = f"C{p_axis[1]}" if p_axis else "None"
                        
            # This calls whatever name you gave your point group classifier method
            group_label = mol.FindPointGroup() 
            
            print(f"{formula:<10} | {planes:<6} | {str(inv):<5} | {p_str:<10} | {group_label}")
        except Exception as e:
            print(f"Error testing {formula}: {e}")