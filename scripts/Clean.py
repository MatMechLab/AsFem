#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 14:01:26 2018
This script will clean the unnecessary file for release
i.e. .iead, cmake-*, .vtu, pvd, *.o
@author: Y. Bai
"""
import os 
import shutil

currentdir=os.getcwd()
print('We are in folder:%s\n'%(currentdir))

ASFEM=0;cmakefolder=0;ideafolder=0;o=0;vtu=0
metafile=0;csvfile=0;valgrind=0;swp=0
cmake=0;cache=0
for subdir,dirs,files in os.walk(currentdir):
    if ('external/eigen' in subdir) or ('external/nlohmann' in subdir) or ('.git' in subdir) or ('.github' in subdir) or ('figures' in subdir):
        continue
    #>>> clean files
    for file in files:
        if ('.csv' in file):
            try:
                csvfile+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif '.swp' in file:
            name=str(file)
            if '.swp' in name[-4:]:
                try:
                    swp+=1
                    removepath=subdir+'/'+file
                    os.remove(removepath)
                except:
                    print('%s is not here'%(file))
        elif ('.json' in file) or ('.cpp' in file) or ('.C' in file) or ('.c' in file and 'cmake_install.cmake' not in file and 'CTestTestfile.cmake' not in file) or ('.h' in file) or ('.hpp' in file) or ('.msh' in file) or ('.geo' in file) or ('.gmsh2' in file) or ('.inp' in file) or ('.py' in file) or ('.C' in file) or ('.txt' in file and 'CMakeCache.txt' not in file) or ('.tex' in file) or ('.jpg' in file) or ('.jpeg' in file) or ('.png' in file) or ('.gif' in file) or ('.pdf' in file) or ('.doc' in file) or ('.docx' in file) or ('.f03' in file) or ('.f08' in file) or ('.f90' in file) or ('.f' in file) or ('.xlsx' in file) or ('Doxyfile' in file):
            continue
        elif ('ASFEM' in file) or ('asfem' in file) or ('asfem-test' in file):
            try:
                ASFEM+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif 'vgcore.' in file:
            try:
                valgrind+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif 'compile_commands.json' in file:
            try:
                cmake+=1
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
        elif ('CMakeCache.txt' in file) or ('cmake_install.cmake' in file) or ('Makefile' in file) or ('CTestTestfile.cmake' in file):
            try:
                cmake+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif ('.vtu' in file) or ('.vtk' in file) or ('.pvd' in file):
            try:
                vtu+=1
                removepath=subdir+'/'+file
                os.remove(removepath)
            except:
                print('%s is not here'%(file))
        elif ('.avi' in file) or ('.gif' in file):
            try:
                metafile+=1
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
        elif '.cache' in dir:
            try:
                cache+=1
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                CacheRemove=True
            except:
                if(not CacheRemove):
                    print('%s is not here'%(dir))
        elif '.clangd' in dir:
            try:
                ideafolder+=1
                removepath=subdir+'/'+dir
                print('remove folder: ',dir)
                shutil.rmtree(removepath)
                ClangRemove=True
            except:
                if(not ClangRemove):
                    print('%s is not here'%(dir))
        elif 'document' in dir:
            try:
                cmakefolder+=1
                removepath=subdir+'/'+dir
                print('remove document folder:',dir)
                shutil.rmtree(removepath)
                DocumentRemove=True
            except:
                if(not DocumentRemove):
                    print('%s is not here'%(dir))
        elif ('cmake-build-debug' in dir) or ('build' in dir) or ('CMakeFiles' in dir) or ('Testing' in dir):
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
print('Remove %4d csv files!'%(csvfile))
print('Remove %4d .o files!'%(o))
print('Remove %4d .swp files!'%(swp))
print('Remove %4d valgrind file!'%(valgrind))
print('Remove %4d meta files!'%(metafile))
print('Remove %4d .idea folder!'%(ideafolder))
print('Remove %4d cache files!'%(cache))
print('Remove %4d cmake file!'%(cmake))
print('Remove %4d cmake folder!'%(cmakefolder))
