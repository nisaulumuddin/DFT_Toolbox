#!/usr/bin/env python

import sys
import os
import numpy as np
import fnmatch 
import math
import argparse, textwrap

sys.path.append('/home/wz300646/Scripts')
import vasp_class_nu as vc


# Upon the need to optimize lattice constants of crystal structures from a website (materials project, etc.), usually the angle components will be correct  but the a, b and c lengths may be different. This script allows the constrained optimization of lattices

# Source of formulas: http://gisaxs.com/index.php/Unit_cell
# In the case of orthorombic structure, where a, b, and c lengths are different..

# MODE:
# --1dac :  Constrain volume, vary c/a ratios (vary c) 
# --1dab : Constrain volume, vary b/a ratios (vary b)
# --vol    : Constrain shape, relax volume (vary multiplier)
# --2dab     : Do a 2D scan of a and b 
# --2dac     : Do a 2D scan of c and a 
# --z      : constrain a and b, vary only along z
# --3d     : do a 3D scan of a, b and c
# --a      : constrain b and c, vary only along a 
# --twinac : 1D scan where a and c are both varied the same value

# SETTING: SCALE/MAG:
# Scale means you have inputted a multiplier to the lattice vectors in current poscar. 
# Mag means you have specified the magnitudes of lengths that will replace the lattice in the current poscar


# OTHER DEGREES OF FREEDOM CAN BE SET UP IN THE INCAR:
# ISIF = 2  : relaxation of atom positions only
# ISIF = 3  : relxation of cell shape, size and atom positions
# ISIF = 4  : relaxation of cell shape and atom positions
# check the VASP manual


# Static calculation version 
# arguments:   [POSCAR]  [MODE]   [A_INI_MAG/SCALE]  [A_FIN_MAG/SCALE]  [STEP_SIZE] [SCALE/MAG] 
# NEED TO FIX THE FORMAT OF THE USER INPUT... NEED TO SPECIFY MODE FIRST THEN THE USER WILL PUT THE CORRECT ARGUMENTS   


parser= argparse.ArgumentParser(description='Optimization of lattice constant through different schemes')
parser.add_argument('-m','--mode', type=str, help= textwrap.dedent('''(REQUIRED) 1dac, 1dab, vol,a,  2dab , 2dac,3d, c, twinac ''') , metavar='', required=True )
parser.add_argument('-i', '--ini', type=float, metavar='', required=False, help='Initial value of scanned parameter for 1dac (c), 1dab (b) and vol (multiplier) modes ')
parser.add_argument('-f','--fin', type=float, metavar='', required=False, help='Final value of scanned parameter for 1dac (c) , 1dab (b), and vol (multiplier) modes')
parser.add_argument('-s','--step', type=int , metavar='', required=True, help='(REQUIRED) Number of steps between initial and final values')
parser.add_argument('-p','--poscar', type=str,metavar='', required=True, help='(REQUIRED) POSCAR file')
parser.add_argument('-t', '--typ', type=str , metavar='', required=False, help='Scaled (0-1) vs magnitude: the type of ini and fin values for 1dac and 1dab modes')


args = parser.parse_args()

def main():
    # ALL THE REQUIRED ARGUMENTS 
    mode = args.mode
    steps = args.step
    poscar = vc.POSCAR(args.poscar)
    


    match mode:
        case "1dac": # a is varied, c is unchanged
            print("a is varied, c is unchanged")

            # ADDITIONAL REQUIRED ARGS 
            ini = args.ini
            fin = args.fin
            typ = args.typ
            
            if typ.startswith("s") or typ.startswith("S"):
                dof_ini = float(ini)*poscar.a
                dof_fin = float(fin)*poscar.a
            elif typ.startswith("m") or typ.startswith("M"):
                dof_ini = float(ini)
                dof_fin = float(fin)
            else:
                print("incorrect input for typ  tag. scale/mag")

            dof_range = np.linspace(dof_ini, dof_fin, steps)
            case_acvol(poscar,dof_range)

        case "1dab": # a is varied, b is unchanged
            print("a is varied, b is unchanged")
            # ADDITIONAL REQUIRED ARGS 
            ini = args.ini
            fin = args.fin
            typ = args.typ

            if typ.startswith("s") or typ.startswith("S"):
                dof_ini = float(ini)*poscar.a
                dof_fin = float(fin)*poscar.a
            elif typ.startswith("m") or typ.startswith("M"):
                dof_ini = float(ini)
                dof_fin = float(fin)
            else:
                print("incorrect input for type tag. scale/mag")
            dof_range = np.linspace(dof_ini, dof_fin, steps)
            case_abvol(poscar,dof_range)

        case "vol": # volume multiplier varied
            print("Volume multiplier is varied")
            # ADDITIONAL REQUIRED ARGS
            ini = args.ini
            fin = args.fin
            
            dof_range = np.linspace(ini,fin,steps) 
            case_vol(poscar, dof_range)

        case "2dab": # Do a 2D scan of a and b
            print("2D scan along a and b")
            print("Additional argument inputs are requested: ")
            case_ab(poscar,steps)
        case "2dac": # Do a 2D scan of c and a
            print("2D scan along c and a")
            print("Additional argument inputs are requested: ")
            case_ac(poscar,steps)

        case "c":
            print("1D scan along c")
            case_c(poscar,steps)

        case "a":
            print("1D scan along a")
            case_a(poscar,steps)

        case "twinac":
            print("1D scan along twin a and c: ")
            case_twinac(poscar,steps)

        
        case "3d":
            print("3D scan along a, b, c")
            print("Additional argument inputs are requested: ")
            case_3d(poscar,steps)
# ----------------------------------------------------------------------



## FUNCTIONS

def case_3d(poscar,step):
    amin = float(input('Insert lower a limit in Angstrom: '))
    amax = float(input('Insert upper a limit in Angstrom: '))
    
    bmin = float(input('Insert lower b limit in Angstrom: '))
    bmax = float(input('Insert upper b limit in Angstrom: '))
    
    cmin = float(input('Insert lower c limit in Angstrom: '))
    cmax = float(input('Insert upper c limit in Angstrom: '))
    
    poscar.lattice_scalar = 1.0
    gamma = poscar.gamma
    alpha = poscar.alpha
    beta = poscar.beta
    
    for a in np.linspace(amin,amax,step):
        for c in np.linspace(cmin,cmax,step):
            for b in np.linspace(bmin,bmax,step):
                factored_lattice_parameters= []
                a_vec = [a, 0.0, 0.0]
                b_vec = [b*np.cos(gamma),b*np.sin(gamma), 0.0]
                c_vec = [c*np.cos(beta),  c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), c*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]
                factored_lattice_parameters = np.array([a_vec, b_vec, c_vec])
                poscar.factored_lattice_parameters  = factored_lattice_parameters
                filename = "POSCAR_a" + str("{:.3f}".format(a)) + "_b" + str("{:.3f}".format(b)) + "_c" + str("{:.3f}".format(c))
                poscar.write_out(filename)

def case_a(poscar,step):
    cmin = float(input('Insert lower a limit in Angstrom: '))
    cmax = float(input('Insert upper a limit in Angstrom: '))

    poscar.lattice_scalar = 1.0
    gamma = poscar.gamma
    alpha = poscar.alpha
    beta = poscar.beta

    a = poscar.a
    b = poscar.b
    c = poscar.c
    b_vec = [b*np.cos(gamma),b*np.sin(gamma), 0.0]
    c_vec = [c*np.cos(beta),  c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), c*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]

    for a in np.linspace(cmin,cmax,step):
        a_vec = [a, 0.0 , 0.0]
        factored_lattice_parameters = np.array([a_vec, b_vec, c_vec])
        poscar.factored_lattice_parameters  = factored_lattice_parameters
        filename = "POSCAR_a" + str("{:.3f}".format(a))
        poscar.write_out(filename)

def case_twinac(poscar,step):
    cmin = float(input('Insert lower a and c limit in Angstrom: '))
    cmax = float(input('Insert upper a and c limit in Angstrom: '))

    poscar.lattice_scalar = 1.0
    gamma = poscar.gamma
    alpha = poscar.alpha
    beta = poscar.beta

    a = poscar.a
    b = poscar.b
    c = poscar.c
    b_vec = [b*np.cos(gamma),b*np.sin(gamma), 0.0]

    for a in np.linspace(cmin,cmax,step):
        a_vec = [a, 0.0 , 0.0]
        c_vec = [a*np.cos(beta),  a*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), a*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]
        factored_lattice_parameters = np.array([a_vec, b_vec, c_vec])
        poscar.factored_lattice_parameters  = factored_lattice_parameters
        filename = "POSCAR_ac" + str("{:.3f}".format(a))
        poscar.write_out(filename)



def case_c(poscar,step):
    cmin = float(input('Insert lower c limit in Angstrom: '))
    cmax = float(input('Insert upper c limit in Angstrom: '))
    
    poscar.lattice_scalar = 1.0
    gamma = poscar.gamma
    alpha = poscar.alpha
    beta = poscar.beta
    
    a = poscar.a
    b = poscar.b
    a_vec = [a, 0.0, 0.0]
    b_vec = [b*np.cos(gamma),b*np.sin(gamma), 0.0]

    for c in np.linspace(cmin,cmax,step):
        c_vec = [c*np.cos(beta),  c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), c*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]
        factored_lattice_parameters = np.array([a_vec, b_vec, c_vec])
        poscar.factored_lattice_parameters  = factored_lattice_parameters
        filename = "POSCAR_c" + str("{:.3f}".format(c))
        poscar.write_out(filename)


def case_vol(poscar,dof_range): #change only the multiplier
    for mult in dof_range:
        poscar.lattice_scalar = mult
        newfilename = 'POSCAR_'+ str(vol) + '_' + str(mult)
        poscar.write_out(newfilename)


def case_ac(poscar,step): # 2D scan of c and a

    amin = float(input('Insert lower a limit in Angstrom: '))
    amax = float(input('Insert upper a limit in Angstrom: '))
    cmin = float(input('Insert lower c limit in Angstrom: '))
    cmax = float(input('Insert upper c limit in Angstrom: '))


    poscar.lattice_scalar = 1.0
    gamma = poscar.gamma
    alpha = poscar.alpha
    beta = poscar.beta 
    b = poscar.b

    constrain = str(input('Is b constrained to equal a or c? (ans: c/a/no ): '))

    for a in np.linspace(amin,amax,step):
        for c in np.linspace(cmin,cmax,step):
            factored_lattice_parameters= []
            a_vec = [a, 0.0, 0.0]
            if constrain == 'a':
                b_vec = [a*np.cos(gamma),a*np.sin(gamma), 0.0]
            elif constrain == 'c':
                b_vec = [c*np.cos(gamma),c*np.sin(gamma), 0.0]
            else:
                b_vec = [b*np.cos(gamma),b*np.sin(gamma), 0.0]
            c_vec = [c*np.cos(beta),  c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), c*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]
            factored_lattice_parameters = np.array([a_vec, b_vec, c_vec])
            poscar.factored_lattice_parameters  = factored_lattice_parameters
            filename = "POSCAR_a" + str("{:.3f}".format(a)) + "_c" + str("{:.3f}".format(c))
            poscar.write_out(filename)




def case_ab(poscar,step):  #2d scan of a and b

    amin = float(input('Insert lower a limit in Angstrom: '))
    amax = float(input('Insert upper a limit in Angstrom: '))
    bmin = float(input('Insert lower b limit in Angstrom: '))
    bmax = float(input('Inser upper b limit in Angstrom: '))

    poscar.lattice_scalar = 1.0
    gamma = poscar.gamma
    beta = poscar.beta
    alpha = poscar.alpha
    c = poscar.c

    constrain = str(input('Is c constrained to equal a or b? (ans: a/b/no ): '))
    
    for a in np.linspace(amin,amax,step):
        for b in np.linspace(bmin,bmax,step):
            factored_lattice_parameters= []
            a_vec = [a, 0.0, 0.0]
            b_vec = [b*np.cos(gamma),b*np.sin(gamma), 0.0]
            if constrain == 'a':                 
                c_vec = [a*np.cos(beta),  a*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), a*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]
            elif constrain == 'b':
                c_vec = [b*np.cos(beta),  b*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), b*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]
            else:
                c_vec = [c*np.cos(beta),  c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), c*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]


            factored_lattice_parameters = np.array([a_vec, b_vec, c_vec])
            poscar.factored_lattice_parameters  = factored_lattice_parameters
            filename = "POSCAR_a" + str("{:.3f}".format(a)) + "_b" + str("{:.3f}".format(b))
            poscar.write_out(filename)


def case_acvol(poscar, dof_range): # change in c/a ratio
    a = poscar.a
    b = poscar.b
    c = poscar.c
    alpha = poscar.alpha
    beta = poscar.beta
    gamma = poscar.gamma    



    constrain2 = str(input('c will be changing according to you initial and final values. do you want volume to be conserved? (ans: y/n) : '))
    
    constrain = str(input('Is b constrained to equal a or c? (ans: c/a/no ): '))
    
    lat_vec = poscar.factored_lattice_parameters
    for c in dof_range: # c is changing according to initial and final values
        opt_vol = float(poscar.volume)


        if constrain2 == 'y' and constrain != 'a' : # b != a
            new_a = opt_vol/(b*c* np.sqrt( 1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2    )  )
        elif constrain2 =='y' and constrain == 'a': # b = a 
            new_a = np.sqrt(opt_vol/(c*np.sqrt( 1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2    )  ))
        else:
            new_a = a 
        
        if constrain == 'c':
            new_b = c
        elif constrain == 'a':
            new_b = new_a
        else:
            new_b = b


        lat_vec[0] = [new_a, 0, 0 ]
        lat_vec[1] = [new_b*np.cos(gamma), new_b*np.sin(gamma), 0]

        lat_vec[2] = [c*np.cos(beta),  c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), c*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]

        ca_ratio = c/new_a

        new_vol = new_a*new_b*c*np.sqrt( 1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2    )


        poscar.lattice_scalar = 1.0
        poscar.factored_lattice_parameters = lat_vec
        newfilename = 'POSCAR_c'+ str("{:.3f}".format(c)) + '_ca' + str("{:.2f}".format(ca_ratio))
        poscar.write_out(newfilename)

        if np.isclose(new_vol, opt_vol, rtol=1e-05, atol=1e-08, equal_nan=False): 
            print("volume is conserved")
        else:
            print("volume is not conserved")



def case_abvol(poscar, dof_range): # b is changing according to the user-input initial and final value
    
    b = poscar.b
    c = poscar.c
    alpha = poscar.alpha
    beta = poscar.beta
    gamma = poscar.gamma 
    lat_vec = poscar.factored_lattice_parameters
    
    constrain2 = str(input('b will be changing according to you initial and final values. do you want volume to be conserved? (ans: y/n) : '))
    constrain = str(input('Is c constrained to equal a or b? (ans: b/a/no ): '))


    
    for b in dof_range: # b is changing

        opt_vol = float(poscar.volume)



        
        if constrain2 == 'y' and constrain != 'a' : # c != a
            new_a = opt_vol/(b*c* np.sqrt( 1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2    )  )
        elif constrain2 =='y' and constrain == 'a': # c = a 
            new_a = np.sqrt(opt_vol/(b* np.sqrt( 1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2    )  ))
        else:
            new_a = a 

        if constrain == 'b':               
            new_c = b
        elif constrain =='a':
            new_c = new_a
        else:
            new_c = c 


        
        lat_vec[0] = [new_a, 0, 0 ] 
        lat_vec[1] = [b*np.cos(gamma), b*np.sin(gamma), 0] #  b is changing 
        lat_vec[2] = [new_c*np.cos(beta),  new_c*((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)), new_c*np.sqrt(1-(np.cos(beta)**2)-(( (np.cos(alpha) - np.cos(beta)*np.cos(gamma))/np.sin(gamma)   )**2))]
        ba_ratio = b/new_a

        new_vol = new_a*b*new_c*np.sqrt( 1 + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2    )


        poscar.lattice_scalar = 1.0
        poscar.factored_lattice_parameters = lat_vec
        newfilename = 'POSCAR_b'+ str("{:.3f}".format(b)) + '_ba' + str("{:.2f}".format(ba_ratio))
        poscar.write_out(newfilename)

        if np.isclose(new_vol, opt_vol, rtol=1e-05, atol=1e-08, equal_nan=False):
            print("volume is conserved")
        else:
            print("volume is not conserved")

if __name__ == '__main__':
        main()

