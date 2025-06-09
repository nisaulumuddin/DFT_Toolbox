#!/usr/bin/env python3

import sys
import os
import numpy as np
import fnmatch 
import math

from vaspreader import vasp_class_nu as vc

# This script is aimed to reorder atomic positions in descending order
# how to use:
# script [POSCAR]  

def main():
        print("")
        poscar = vc.POSCAR(sys.argv[1])
        poscar.reorder_coord('d') # order descending z-coord
        poscar.write_out(f'POSCAR_reordered')
if __name__ == '__main__':
        main()

