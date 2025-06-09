#!/usr/bin/env python3

import numpy as np
import sys

sys.path.append('/home/wz300646/Scripts/Utilities_nu/')
import vasp_class_nu as vc



incar = vc.INCAR('INCAR')

incar.set_tag('ISIF','2')
incar.set_tag('EDIFF','1E-6')
incar.set_tag('EDIFFG','-0.02')
incar.set_tag('LREAL','False')
incar.set_tag('ADDGRID','True')
incar.set_tag('ISYM','0')
incar.set_tag('NSW','0')
incar.set_tag('IBRION','-1')
incar.set_tag('PREC','Accurate')

incar.write_out('INCAR_s2c')

