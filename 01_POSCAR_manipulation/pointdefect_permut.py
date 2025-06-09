#!/usr/bin/env python3



# This script permutes the coordinates of selected atoms. good for site occupancy exploration given a fixed composition

from itertools import combinations, tee
import datetime
import numpy as np
import sys
import os
import itertools as it
from vaspreader import vasp_class_nu as vc
import operator
import copy
import shutil
from itertools import combinations
import operator
from concurrent.futures import ProcessPoolExecutor, as_completed


from pymatgen.core.structure import Structure
from pymatgen.analysis.structure_matcher import StructureMatcher

def lexico_permute(a):
    a = list(a)
    yield a
    n = len(a) - 1
    while True:
        for j in range(n-1, -1, -1):
            if a[j] < a[j + 1]:
                break
        else:
            return

        v = a[j]
        for k in range(n, j, -1):
            if v < a[k]:
                break

        a[j], a[k] = a[k], a[j]
        a[j+1:] = a[j+1:][::-1]
        yield a

def check_structure(i, permut_i, selected_pos, unchanged_index, unchanged_atoms,
                    section_index, poscar_data, lattice_parameters, atom_list, structure0_data):
    from pymatgen.core.structure import Structure
    from pymatgen.analysis.structure_matcher import StructureMatcher
	

    matcher = StructureMatcher(primitive_cell=False)  # ✅ move inside

    # Recreate structure0
    structure0 = Structure(lattice=lattice_parameters, species=atom_list, coords=structure0_data)

    # Prepare current structure
    zipped = list(zip(permut_i, selected_pos))
    sorted_zipped = sorted(zipped, key=operator.itemgetter(0))
    _, sorted_list2 = zip(*sorted_zipped)

    fin_pos = np.zeros_like(poscar_data)
    fin_pos[unchanged_index] = unchanged_atoms
    fin_pos[section_index] = sorted_list2

    structure = Structure(lattice=lattice_parameters, species=atom_list, coords=fin_pos)

    # Compare structure
    if matcher.fit(structure0, structure):
        return None  # Duplicate
    else:
        return structure  # Return non-duplicate for further checking



if __name__ == "__main__":
    
    print("Always check your POSCAR structures before running. Adjust order of POTCAR accordingly.")
    poscar =  vc.POSCAR(sys.argv[1]) # load the file
    print(poscar.atom_header)
    user = str(input("write the list of elements that you wish to permute separated by space i.e. Fe Al:     "))
    print( datetime.datetime.now())
    print(user)

    # find the coordinates of the user-specified elements in the POSCAR 

    len_permut = 0
    len_unchanged=0
    atom_list = []
    for element, quantity in zip(poscar.atom_header[0],poscar.atom_header[1]):
        if element in user:
            len_permut+=quantity
            atom_list.extend([element]*quantity)
        else:
            len_unchanged += quantity

    order = sorted([atom[i:i+2] for atom in atom_list for i in range(0, len(atom), 2)]) # this needs to be in alphabetical order for the lexico_permute function to work

    selected_atoms = np.zeros((len_permut,3))
    unchanged_atoms = np.zeros((len_unchanged,3))
    i=0
    j=0
    index=0
    section_index = []
    unchanged_index = []

    for element, pos_data in zip(poscar.atom_list, poscar.pos_data):
        if element in user:
            selected_atoms[i]= pos_data
            i+=1
            
            section_index.append(index)
            index+=1
        else:
            unchanged_atoms[j] = pos_data
            j+=1
            
            unchanged_index.append(index)
            index+=1

    selected_pos = selected_atoms.copy() # section of POSCAR in which atom coordinates should be permuted
    unchanged = unchanged_atoms.copy()

    # permutation function 


    permut = []
    for i, u in enumerate(lexico_permute(order), 0):
        # print(i,u)
        permut.append(list(u))


    atom_list = []
    for speciesnum,species in enumerate(poscar.atom_header[0]):
            for occur in range(poscar.atom_header[1][speciesnum]):
                    atom_list.append(species)
    poscar.atom_list = atom_list

    user_input = str(input("Do you want equivalent POSCARs to be erased? (y/n): "))

    if user_input.startswith("y") or user_input.startswith("Y"):
        print("putting permuted coordinates into the POSCAR")
        print( datetime.datetime.now())
        list_all_poscar = []

        # putting the permuted coordinates into the POSCAR
        matcher = StructureMatcher(primitive_cell=False)

        # Function to check if a structure is a duplicate of any structure in a list.
        def is_duplicate(new_struct, struct_list, matcher):
            for s in struct_list:
                if matcher.fit(new_struct, s):
                    return True
            return False

        # Outside the function
        unique_structures = []

        structure0_data = copy.deepcopy(poscar.pos_data)  # because structure0 is constant
        poscar_data = copy.deepcopy(poscar.pos_data)

        with ProcessPoolExecutor() as executor:
            futures = []
            for i in range(len(permut)):
                futures.append(executor.submit(
                    check_structure,
                    i, permut[i], selected_pos, unchanged_index, unchanged_atoms,
                    section_index, poscar_data, poscar.lattice_parameters, poscar.atom_list, structure0_data
                ))

            for future in as_completed(futures):
                result = future.result()
                if result is not None:
                    if not is_duplicate(result, unique_structures, matcher):
                        unique_structures.append(result)
                    
            
        index = 0 
        for structure in unique_structures:
            filename = "POSCAR_" + str(index) + ".vasp"
            structure.to(fmt="poscar", filename=filename)
            index+=1
    else:
        for i in range(len(permut)):
            zipped=list(zip(permut[i],selected_pos))
            sorted_zipped=sorted(zipped, key=operator.itemgetter(0))
            sorted_list1, sorted_list2 = zip(*sorted_zipped) # zipped : element , coordinates
            fin_pos = np.zeros((len(poscar.pos_data),3))
            fin_pos[unchanged_index] = unchanged_atoms
            fin_pos[section_index] = sorted_list2
            poscar.pos_data = copy.deepcopy(fin_pos)
            filename = "POSCAR_" + str(i) + ".vasp"
            poscar.write_out(filename)
        
    print( datetime.datetime.now()) 
