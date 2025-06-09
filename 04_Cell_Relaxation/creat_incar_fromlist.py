#!/usr/bin/env python3

import numpy as np
import sys
import argparse, textwrap

from vaspreader import vasp_class_nu as vc



def main(incar,text,mode):

    #args = parser.parse_args() 


    Incar = vc.INCAR(incar)
    inpfile = open(text,'r')
    inpcont = inpfile.readlines()

    count = 0

    if mode == "new":
        run = uniq_new(inpcont,Incar)

    elif mode == "all":
        run = insrt_all(inpcont,Incar)
    else:
        print("invalid mode")

def insrt_all(inpcont,Incar):
    for line in inpcont:
        tag_info = line.split()
        tag_name = tag_info[0]
        tag_det = " ".join(str(element) for element in tag_info[2:])
        Incar.set_tag(tag_name, tag_det)

    filename = "INCAR_insrt" 
    Incar.write_out(filename)


def  uniq_new(inpcont,Incar):
    count = 0
    for line in inpcont:
        tag_info = line.split()
        tag_name = tag_info[0]
        tag_det = " ".join(str(element) for element in tag_info[2:])
        Incar.set_tag(tag_name, tag_det)
        filename = "INCAR_" + str(count)
        Incar.write_out(filename)
        count+=1


if __name__ == "__main__":

    parser= argparse.ArgumentParser(description='Replaces/Adds a line the INCAR file according to a text file (new INCAR made for every line in the text file)')
    parser.add_argument('-i', '--incar', type=str, help= 'INCAR file' , metavar='', required=True )
    parser.add_argument('-t', '--text', type=str, metavar='',  help='text file with list of INCAR tags to replace with from the original INCAR', required=True)
    parser.add_argument('-m', '--mode', type=str, help= 'mode: new (new INCAR made for every line in the text file) vs all (all lines in text file are implemented into one new INCAR)' , metavar='', required=True )
    
    args = parser.parse_args()
    main(incar=args.incar, text=args.text, mode =args.mode)
