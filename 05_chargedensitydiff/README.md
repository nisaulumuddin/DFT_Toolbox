# How to run a differential charge density analysis for a bulk crystal

## Step 0 : Reoptimize structure with high k-points using the following setup 
```
LREAL = FALSE
LORBIT = 10
ISMEAR = -5
IBRION = -1
NSW = 0
LCHARG = TRUE 
LAECHG = TRUE
```


## Step 1: Run script
```
usage: setup_chgdiff.py [-h] --targetdir TARGETDIR --inidir INIDIR

Your script description here

options:
  -h, --help            show this help message and exit
  --targetdir TARGETDIR
                        The target directory for the operation
  --inidir INIDIR       The directory to extract charge density from (the directory in which step1 was done)


setup_chgdiff.py --targetdir <target_dir> --inidir <ini_dir>
```


## Step 2: Run the calculations in the cluster and extract charge density difference in VESTA 

