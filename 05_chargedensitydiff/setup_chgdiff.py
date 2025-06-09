#!/usr/bin/env python3

import numpy as np
import sys

from vaspreader import vasp_class_nu as vc
import vasp6_class_nu as vc6
import argparse
import os
import glob

import shutil


# Any custom imports
# from mymodule import my_function

def main(args):
    """
    Provide the directory from which charge density will be calculated and a target directory for the calculation jobs.
    Run calculations and use the CHGCAR to calculate charge density difference
    """
    ini_dir = args.inidir
    target_dir = args.targetdir

    poscar_path = os.path.join(ini_dir, 'POSCAR')
    incar_path = os.path.join(ini_dir, 'INCAR')
    outcar_path = os.path.join(ini_dir, 'OUTCAR')
    
    poscar = vc6.POSCAR(poscar_path)
    incar = vc6.INCAR(incar_path)
    outcar = vc6.OUTCAR(outcar_path)

    potcar_i = []
    for element in poscar.atom_header[0]:
        potcar_i.append(input(f"Provide the path of the POTCAR file of {element}: "))
    potcar_path = list(zip(poscar.atom_header[0],potcar_i))

    # creating POSCAR files, POTCAR files, and copy KPOINT files
    setup_poscar(poscar)
    setup_potcar(potcar_path, target_dir)
    copy_kpoints(ini_dir, target_dir)

    poscar_files = [
        f for f in os.listdir('.')
        if f.startswith('POSCAR_') and f.endswith('.vasp') and os.path.isfile(f)
    ]

    filenames = [filename.removeprefix('POSCAR_').removesuffix('.vasp') for filename in poscar_files]

    # creating the INCAR files
    for filename, poscar_i in zip(filenames,poscar_files):
        poscar_i = vc.POSCAR(poscar_i)
        edit_incar(incar, poscar_i, outcar, filename)

    # move INCAR and POSCAR files 
    mv_incar_poscar(target_dir)


def parse_args():
    """
    Parse command-line arguments using argparse.
    """
    parser = argparse.ArgumentParser(description="Your script description here")

    # Define the arguments you want to accept
    parser.add_argument(
        '--targetdir', 
        type=str, 
        required=True, 
        help="The target directory for the operation"
    )

    parser.add_argument(
        '--inidir', 
        type=str, 
        required=True, 
        help="The directory to extract charge density from"
    )


    return parser.parse_args()

def setup_poscar(poscar):
    zip_atoms = zipped = zip(poscar.atom_list , poscar.pos_data)

    count=0
    for atom,pos in zip_atoms:
        # create POSCAR
        poscar.atom_header = [[atom],[1]]
        poscar.pos_data = pos
        poscar.selective_dynamics= False
        poscar.fixation_data = None

        newname = 'POSCAR_' + str(atom) + "_" + str(count) + ".vasp"
        poscar.write_out(newname)
        count+=1

def setup_potcar(potcar_path, parent_dir):
        # create folders
        # Remove parent folder if it exists and recreate
        if os.path.exists(parent_dir): # the chgdif folder
            shutil.rmtree(parent_dir) # the chgdif folder
        os.makedirs(parent_dir)

        poscar_files =  glob.glob('POSCAR_*')
        for file in poscar_files: 
            identifier = file.removeprefix('POSCAR_').removesuffix('.vasp')
            atom = identifier.split('_')[0]
            # Define target directory
            target_dir = os.path.join(parent_dir, identifier)
            os.makedirs(target_dir, exist_ok=True)

            # copy the right POTCAR 
            for element, potcar_path_i in potcar_path:
                if element == atom:
                    # Check if it's a directory
                    if os.path.isdir(target_dir):
                        # Copy KPOINTS into the folder
                        shutil.copy(potcar_path_i, os.path.join(target_dir, 'POTCAR'))
                        break


def mv_incar_poscar(parent_dir):

    # Find all INCAR files matching the pattern
    incar_files = glob.glob('INCAR_*')

    for incar in incar_files:
        # Extract the identifier (e.g., Ta_1)
        identifier = incar.replace('INCAR_', '')
        # Build corresponding POSCAR filename
        poscar = f'POSCAR_{identifier}' + ".vasp"
        # Define target directory
        target_dir = os.path.join(parent_dir, identifier)
        # Copy both files into the target directory
        shutil.move(incar, os.path.join(target_dir, 'INCAR'))
        if os.path.exists(poscar):
            shutil.move(poscar, os.path.join(target_dir, 'POSCAR'))
        else:
            print(f'Warning: {poscar} not found.')
    print('INCAR and POSCAR files are in separate folders!')

def copy_kpoints(inifolder, parent_dir):
    print("Copying KPOINTS to each target folder ")
    # Path to the source KPOINTS file
    kpoints_file = os.path.join(inifolder, 'KPOINTS')

    # Check if KPOINTS file exists
    if not os.path.exists(kpoints_file):
        print(f'Error: {kpoints_file} not found.')
    else:
        # Loop over subfolders inside target parent folder
        for folder_name in os.listdir(parent_dir):
            folder_path = os.path.join(parent_dir, folder_name)

            # Check if it's a directory
            if os.path.isdir(folder_path):
                # Copy KPOINTS into the folder
                shutil.copy(kpoints_file, folder_path)
    print(f'Copied KPOINTS to each target folder')


def edit_incar(incar, poscar, outcar, filename):
    atom_header = poscar.atom_header
    elements, counts = atom_header
    rwigs_list = list(outcar.element_rwigs)
    rwigs_dict = dict(rwigs_list)

    rwigs_string = []
    stradd = []
    atom_list = []
    # edit this yourself
    for element, count in zip(elements,counts):
        if element == 'Ta':
            Fe_multiplier = count
            Fe_mag = str('0 0 0 ')
            str_add = int(Fe_multiplier)*str(Fe_mag)
            stradd.append(str(str_add))
            rwigs_string.append(rwigs_dict['Ta'])
            atom_list.append(element)
        elif element == 'Fe':
            Ta_multiplier = count
            Ta_mag = str('0 0 1.5 ')
            str_add = int(Ta_multiplier)*str(Ta_mag)
            stradd.append(str(str_add))
            rwigs_string.append(rwigs_dict['Fe'])
            atom_list.append(element)

    m_constr = ' '.join(stradd)

    rwigs_string = ' '.join(rwigs_string)
    
    incar.set_tag('ISMEAR','-5')
    incar.set_tag('EDIFF','1E-6')
    incar.set_tag('LREAL','FALSE')
    incar.set_tag('LCHARG','TRUE')
    incar.set_tag('LORBIT','10')
    incar.set_tag('LREAL','FALSE')
    incar.set_tag('LAMBDA','1')
    incar.set_tag('RWIGS', rwigs_string)
    incar.set_tag('I_CONSTRAINED_M','1')
    incar.set_tag('ISPIN','2')
    incar.set_tag('M_CONSTR',str(m_constr))
    incar.set_tag('LAECHG','TRUE')
    incar.set_tag('NSW','0')
    incar.set_tag('IBRION','-1')
    incar.write_out('INCAR_' + str(filename))
    return atom_list

if __name__ == '__main__':
    # Parse the arguments
    args = parse_args()

    # Call the main function with the parsed arguments
    main(args)
