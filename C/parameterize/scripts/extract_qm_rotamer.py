#!/usr/bin/env python

import os

os.system('/usr/bin/find ./ -name "output.out" > find-qm-list.txt')

fIn = open('find-qm-list.txt','r')

for line in fIn:
    str.split('/')
    strlist=line.split('/')
    os.system('cp ./'+strlist[1]+'/output.out '+'qm-'+strlist[1]+'.out')
    os.system('cp ./'+strlist[1]+'/mol-qm-rotamer.gjf '+'qm-'+strlist[1]+'.gjf')

fIn.close()

