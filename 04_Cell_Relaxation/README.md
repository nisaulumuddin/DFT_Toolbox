# Workflow for cell relaxation in VASP 

Based on our experience with VASP, the simulation cell (a,b,c lengths) relaxation scheme implemented in VASP does not always converge to structures with minimal external stress. The presence of external stress for the x,y and z components indicate that an external pressure is required to keep the structure in equilibrium despite it being a local minimum structure. 

To do simulation cell relaxation, it is better to use VASP 6.4.3. This will enable additional constraint modes allowed in VASP. 

## Step 1: Use VASP algorithm to converge to a local minimum structure

Do a few cycles of the (cell + atomic position) relaxation scheme followed by the (atomic position) relaxation scheme.

(cell + atomic position) relaxation
```
PREC = Accurate # needs to be high enough
ENCUT = 550 # needs to be high enough for the system
ISIF = 8 
EDIFF = 1E-6
POTIM = 0.2
NSW = 50
LATTICE_CONSTRAINTS = .TRUE. .TRUE. .TRUE. # can adjust if you want to relax only in certain directions
```
(atomic position) relaxation
```
PREC = Accurate # needs to be high enough
ENCUT = 550 # needs to be high enough for the system
ISIF = 2
EDIFF = 1E-6
POTIM = 0.1
NSW = 100 
LATTICE_CONSTRAINTS = .TRUE. .TRUE. .TRUE. # can adjust if you want to relax only in certain directions
```

When cell relaxation occurs for magnetic material systems, the starting magnetic moment can change as you switch from the (cell + atomic position) relaxation to the (atomic position) relaxation scheme. Therefore, it is sometimes necessary to constrain the magnetic configuration either by (direction) or by (direction + magnitude)

A sample of magnetic direction constraint
```
I_CONSTRAINED_M = 1  # direction
M_CONSTR= 3*96*0  0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3 0 0 3
LAMBDA = 1
RWIGS =  1.38 1.25
```
You can use the `incar_tags.py` or the `setup_incar_fromlist.py`script to set up the INCAR tags for multiple directories.

## Step 2: Check the final stresses in the last OUTCAR

You can do this easily by typing

```
grep -B 3 external OUTCAR
```

If the stresses are less than 0.5 kBar, then the structure is well converged. (Maximum tolerable is 1 kBar)
If not, go to Step 3.

## Step 3: Manual scan along the lattice parameter

Use the `latstruct_gen.py` script to run a manual scan over your degrees of freedom specific to your material system

```
usage: latstruct_gen.py [-h] -m  [-i] [-f] -s  -p  [-t]

Optimization of lattice constant through different schemes

options:
  -h, --help      show this help message and exit
  -m , --mode     (REQUIRED) 1dac, 1dab, vol, 2dab , 2dac, c
  -i , --ini      Initial value of scanned parameter for 1dac (c), 1dab (b) and vol (multiplier) modes
  -f , --fin      Final value of scanned parameter for 1dac (c) , 1dab (b), and vol (multiplier) modes
  -s , --step     (REQUIRED) Number of steps between initial and final values
  -p , --poscar   (REQUIRED) POSCAR file
  -t , --typ      Scaled (0-1) vs magnitude: the type of ini and fin values for 1dac and 1dab modes
```
