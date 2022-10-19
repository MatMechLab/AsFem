#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat 2022.01.01
This script will update the link of AsFem homepage
@author: Yang Bai
"""
import os 
import shutil
import fileinput

OldStr='//* https://github.com/yangbai90/AsFem.git'
NewStr='//* https://github.com/M3Group/AsFem'

currentdir=os.getcwd()
print('We are in folder:%s\n'%(currentdir))

nFiles=0;nLines=0
for subdir,dirs,files in os.walk(currentdir):
    if ('external/eigen' in subdir) or ('.git' in subdir):
        continue
    for file in files:
        if (('.h' in file) or ('.cpp' in file) or ('.hpp' in file) or ('.cxx' in file)) and ('.bak' not in file):
            # here we only edit the source code file
            oldfile=subdir+'/'+file
            newfile=subdir+'/'+file+'.bak'
            fin=open(oldfile,'r',encoding="ISO-8859-1") # open old file
            fout=open(subdir+'/'+file+'.bak','w+') # write to old.bak file
            for line in fin:
                fout.write(line.replace(OldStr,NewStr))
                nLines+=1
            fin.close()
            fout.close()
            os.remove(oldfile)
            os.rename(newfile,oldfile)
            nFiles+=1

print('Update all files, %d files have been modified!'%(nFiles))
print('Total lines of AsFem (include comments) code is:%d'%(nLines))
