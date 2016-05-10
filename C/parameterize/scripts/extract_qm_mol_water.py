#!/usr/bin/env python

import os


os.system('/usr/bin/find ./ -name "opt-output.out" > find-qm-list.txt')

fIn = open('find-qm-list.txt','r')

for line in fIn:
    str.split('/')
    strlist=line.split('/')
    newname=strlist[1].replace('dat-','')
    os.system('cp ./'+strlist[1]+'/QM/opt-output.out '+'qm-mol-wat-'+newname+'.out')
    os.system('cp ./'+strlist[1]+'/QM/dimer.gjf '+'qm-mol-wat-'+newname+'.gjf')
    os.system('cp ./'+strlist[1]+'/mol-waters.pdb '+'qm-mol-wat-'+newname+'.pdb')

fIn.close()

