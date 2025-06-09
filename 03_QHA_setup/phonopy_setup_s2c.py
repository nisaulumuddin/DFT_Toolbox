#!/usr/bin/zsh

for f in Fe*
do 
cd $f
cp POSCAR-001 RECELL
/home/wz300646/Installations/AELAS_v1_0/bin/AELAS -ka
cp NEWKPT KPOINTS
    for i in POSCAR-*
    do
    mkdir ${i/POSCAR-/}
    cp $i ${i/POSCAR-/}/POSCAR
    cp INCAR KPOINTS POTCAR ${i/POSCAR-/}
    done
cd .. 
done 
