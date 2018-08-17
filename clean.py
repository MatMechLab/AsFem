#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:01:26 2018
Try to write a simple msh file for MOOSE
@author: Y. Bai
"""
import os 
import shutil

currentdir=os.getcwd()
print('We are in folder:%s\n'%(currentdir))

vtu=0;jit=0;pycache=0
for subdir,dirs,files in os.walk(currentdir):
    #>>> clean files
    for file in files:
        if 'vtu' in file:
            try:
                vtu+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
    #>> clean folder
    for dir in dirs:
        if '.jit' in dir:
            try:
                jit+=1
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                IdeaRemove=True
            except:
                if(not IdeaRemove):
                    print('%s is not here'%(dir))
        if '__pycache__' in dir:
            try:
                pycache+=1
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                IdeaRemove=True
            except:
                if(not IdeaRemove):
                    print('%s is not here'%(dir))


################################################
print('Remove %4d .vtu files!'%(vtu))
print('Remove %4d .jit folder!'%(jit))
print('Remove %4d pycache folder!'%(pycache))
