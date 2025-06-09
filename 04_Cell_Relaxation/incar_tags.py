#!/usr/bin/env python3

import numpy as np
import sys

from vaspreader import vasp_class_nu as vc



#---- extract from POSCAR the number of Fe atoms 

poscar = vc.POSCAR('POSCAR')

atom_header = poscar.atom_header
elements, counts = atom_header

stradd = []
for element, count in zip(elements,counts):
    if element == 'Fe':
        Fe_multiplier = count
        Fe_mag = str('0 0 3 ')
        str_add = int(Fe_multiplier)*str(Fe_mag)
        stradd.append(str(str_add))
    elif element == 'Ta':
        Ta_multiplier = count
        Ta_mag = str('0 0 0 ')
        str_add = int(Ta_multiplier)*str(Ta_mag)
        stradd.append(str(str_add))

m_constr = ' '.join(stradd)


incar = vc.INCAR('INCAR')

incar.set_tag('ISPIN','2')
incar.set_tag('ISIF','8')
incar.set_tag('NSW','50')
incar.set_tag('POTIM','0.2')
incar.set_tag('EDIFF','1E-6')
incar.set_tag('EDIFFG','-0.02')
incar.set_tag('LATTICE_CONSTRAINTS','.FALSE. .FALSE. .TRUE.')
incar.set_tag('I_CONSTRAINED_M','1')
incar.set_tag('LAMBDA','1')
incar.set_tag('RWIGS','1.25 1.38')
incar.set_tag('M_CONSTR',str(m_constr))
incar.write_out('INCAR_isif8')
