#!/usr/bin/env python3

import numpy as np
import sys, os
# Make sure that the VASP classes is in the parent directory's subdirectory 'Utilities'
from vaspreader import vasp_class_nu as vc
import vasp_class_nu as vc
import func_magmomsetup as func

import argparse
import sys

def main(args):
    print(f"Arguments: {args}")
    
    # the purpose of this script is to adjust the MAGMOM tag of  an INCAR file for a supercell calculation
    # works only in the case if atoms that belong to the same atomic layer (within 10E-3 Angstrom accuracy) possess the same magmom

    poscar_cell = vc.POSCAR(args.unitcell) # the initial cell
    poscar_sup = vc.POSCAR(args.supercell) # the supercell
    incar = vc.INCAR(args.incar) # the INCAR of the initial cell 
    supercell_input = args.size # user input

    small_magmomlist = incar.tags['MAGMOM']
    small_header = poscar_cell.atom_header
    supercell = [int(x) for x in supercell_input.split()]

    
    # first step, reorder the supercell in descending order (edit this so that it has a range of acceptance)
    poscar_sup_desc = poscar_sup.reorder_coord('d')

    # user inputs the layer-by-layer , element-by-element magmom list of the initial cell
    magmom_cell_string = args.magmom # user input
    magmom_cell = func.split_string(magmom_cell_string)
    mapped_elem,new_atomheader = func.map_elements(small_header,magmom_cell)
    magmom_magn = magmom_cell[1]
    #magmom_sup= " ".join([f"{int(value)}*{factor}" for value, factor in zip(supercell_quantity, magmom_magn)])

    duplicates = func.find_duplicates(new_atomheader)
    #result['Fe']
    # print(magmom_magn[duplicates['Fe']])

    supercell_quantity= np.zeros(len(small_magmomlist))
    count=0
    for f in mapped_elem[0]:
        element_quantity_layer = f[1]*supercell[0]*supercell[1] # element quantity multiplied by x and y supercell multiplier
        supercell_quantity[count] = element_quantity_layer
        count+=1
    
    magmom_sup = ""
    count = 0
    for elem in new_atomheader:
        if elem in duplicates:
            magmom = [int(magmom_magn[i]) for i in duplicates[elem]]
            quantity = [int(supercell_quantity[i]) for i in duplicates[elem]]
            magmom_sup += " ".join([f"{int(value)}*{factor} " for value, factor in zip(quantity, magmom)])
                        
        if elem not in duplicates:
            magmom_sup += f" {int(supercell_quantity[count])}*{int(magmom_magn[count])} "*supercell[2]
        count+=1
        
    incar.set_tag('MAGMOM',magmom_sup)
    newname = 'INCAR_NEW'
    incar.write_out(newname)
    
    
if __name__ == "__main__":
    # Example: Simulating arguments
    #sys.argv = ["script.py", "--name", "Alice", "--age", "25", "--verbose"]
    
    # Initialize parser
    parser = argparse.ArgumentParser(description="Argparse with sys.argv example.")
    parser.add_argument("--unitcell", "-uc", type=str, required=True, help="POSCAR name of unit cell")
    parser.add_argument("--supercell", "-sc", type=str,required=True, help="POSCAR name of supercell")
    parser.add_argument("--incar", "-i", type=str,required=True, help="INCAR of unit cell")
    parser.add_argument("--size", "-s", type=str,required=True, help="supercell multipliers : provide in x y z")
    parser.add_argument("--magmom", "-m", type=str,required=True, help="atom*magmom list of unit cell in the order of atomic layers")
    
    # Parse the arguments
    args = parser.parse_args()
    
    main(args)

