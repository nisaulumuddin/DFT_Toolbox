# Setup for Quasi Harmonic Approximation in VASP 

We follow the same procedure as https://www.youtube.com/watch?v=FX7WjL074g4

There are 2 main steps:
-  Step 1: structure relaxation
-  Step 2: Phonon calculations (using CONTCAR from Step 1)



## Step 1 : Structure Relaxation 

Set up a calculation as normal but add these additional tags:
```
ISYM = 1 
LREAL = False
PREC = Accurate
ISIF = 2
```

## Step 2: Phonon Calculations

### Step 2a:
 
Here we will use the supercell method (finite-difference method) to isolate atom-atom interactions when doing phonon calculations 

This is where we will use Phonopy (https://phonopy.github.io/phonopy/install.html)

use the CONTCAR structure from Step 1 

We need to creat a supercell where the lattice constants of the simulation cell are at least 1-1.5nm 

To do this, run command: 
```
phonopy -d --dim = "[supercell multipliers]"
```

example: phonopy -d --dim = "4 4 1"

### Step 2b: 

This will create `POSCAR-00#` which each contains single atomic displacements 

### Step 2c: 

Create calculation folders for each POSCAR-00#, and put each POSCAR in its corresponding calculation folder. 

INCAR tags:
```
LREAL = False
ADDGRID = True
ISYM = 0 
NSW = 0 
IBRION = -1
PREC = Accurate 
```

Change KPOINTS according to the new dimensions (keep KPPRA the same)

Scripts to et up these jobs

#### phonopy_magmomsetup_s2c.py 

This script sets up the MAGMOM tag in the INCAR for the SUPERCELL (in POSCAR-###) based on user-given information regarding the supercell dimensions and the layer-by-layer MAGMOM information of the unit cell. The POSCAR atom header should have the same number of entries of the --magmom list provided by the user, so that the script can associate the atom number to its elements.

```
usage: phonopy_magmomsetup_s2c.py [-h] --unitcell UNITCELL --supercell
                                  SUPERCELL --incar INCAR --size SIZE --magmom
                                  MAGMOM

Argparse with sys.argv example.

optional arguments:
  -h, --help            show this help message and exit
  --unitcell UNITCELL, -uc UNITCELL
                        POSCAR name of unit cell
  --supercell SUPERCELL, -sc SUPERCELL
                        POSCAR name of supercell
  --incar INCAR, -i INCAR
                        INCAR of unit cell
  --size SIZE, -s SIZE  supercell multipliers : provide in x y z
  --magmom MAGMOM, -m MAGMOM
                        atom*magmom list of unit cell in the order of atomic
                        layers

```

Example of use:

```
phonopy_magmomsetup_s2c.py -uc POSCAR -sc POSCAR-001 -i INCAR -s '3 3 2' -m '4*0 6*0 1*-3 1*3' 
```
This script generates `INCAR_NEW`

#### phonopy_incar_s2c.py

This script sets up the rest of the INCAR tag needed to run for step 2c (the INCAR tags listed above)

```
Sets up INCAR file for phonopy (step 2c)

optional arguments:
  -h, --help     show this help message and exit
  -i , --incar   INCAR file

```

Example of use:
```
phonopy_incar_s2c.py -i INCAR
```
This script generates `INCAR_s2c` 



### Step 2d: Create Force sets
After running these calculations, extract vasprun.xml from each folder

```
phonopy -f dis-001/vasprun.xml dis-002/vasprun.xml dis-003/vasprun.xml dis-004/vasprun.xml ...
```

A useful bash command, if the folders start with "0":
```
for f in 0*; do echo -n "$f/vasprun.xml " >> dirlist ; done
phonopy -f $(cat dirlist) 
```


## Step 3: Post-processing 
Force constants are calculated from the sets of forces. A part of the dynamical matrix is built from the force constants. Phonon frequences and eigenvectors are calculated from the dynamical matrices with the specified q-points. See [QHA Theory](./phonons_theory.pdf) to read its theory in more depth. 

#### Phonon DOS
Create `mesh.conf`
Example of `mesh.conf`
```
ATOM_NAME = Cr Al C 
DIM = 4 4 1 # the supercell dimensions 
MP = 24 24 48  # the grid points in x y and z dir (get from NGXF, NGYF, NGZF)
```
Then run command:
```
phonopy -p -s mesh.conf
```


#### Thermal Properties 

Run command:
```
phonopy -p -s -t mesh.conf
```
The red like in Hemholtz Free Energy A(T)
```
A(T) = U(T)-TS(T)
A(T) = [Lattice internal energy] + [Vibrational internal energy] - T[Vibrational entropy]
```

#### Phonon Band Structure

Band.conf example
```
ATOM_NAME = Cr Al C
DIM = 4 4 1 #supercell
PRIMITIVE_AXIS = AUTO
BAND = [find the suitable band path for our specific structure]
BAND_POINTS = 101 # how many data points 
BAND_CONNECTION = .TRUE. # line graph
```

```
phonopy -p -s  band.conf
```

Check that there are no imaginary modes (negative) to make sure that your structure is stable dynamically. 
Check https://en.wikipedia.org/wiki/Brillouin_zone for the high-symmetry band path