#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:01:26 2018
This script will clean the unnecessary file for release
i.e. .iead, cmake-*, .vtu, *.o
@author: Y. Bai
"""
import os 
import shutil

currentdir=os.getcwd()
print('We are in folder:%s\n'%(currentdir))

ASFEM=0;cmakefolder=0;ideafolder=0;o=0;vtu=0
cmake=0
for subdir,dirs,files in os.walk(currentdir):
    #>>> clean files
    for file in files:
        if ('.i' in file) or ('.cpp' in file) or ('.C' in file) or ('.c' in file and 'cmake_install.cmake' not in file) or ('.h' in file) or ('.hpp' in file) or ('.msh' in file) or ('.geo' in file) or ('.gmsh2' in file) or ('.py' in file) or ('.C' in file) or ('.txt' in file and 'CMakeCache.txt' not in file) or ('.tex' in file) or ('.jpg' in file) or ('.png' in file) or ('.pdf' in file):
            continue
        elif ('ASFEM' in file) or ('asfem' in file):
            try:
                ASFEM+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif '.o' in file:
            try:
                o+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif ('CMakeCache.txt' in file) or ('cmake_install.cmake' in file) or ('Makefile' in file):
            try:
                cmake+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif ('.vtu' in file) or ('.vtk' in file):
            try:
                vtu+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
    #>> clean folder
    for dir in dirs:
        if '.idea' in dir:
            try:
                ideafolder+=1
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                IdeaRemove=True
            except:
                if(not IdeaRemove):
                    print('%s is not here'%(dir))
        elif ('cmake-build-debug' in dir) or ('build' in dir) or ('CMakeFiles' in dir):
            try:
                cmakefolder+=1
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                IdeaRemove=True
            except:
                if(not IdeaRemove):
                    print('%s is not here'%(dir))

print('Remove %4d ASFEM files!'%(ASFEM))
print('Remove %4d vtu files!'%(vtu))
print('Remove %4d .o files!'%(o))
print('Remove %4d .idea folder!'%(ideafolder))
print('Remove %4d cmake file!'%(cmake))
print('Remove %4d cmake folder!'%(cmakefolder))
