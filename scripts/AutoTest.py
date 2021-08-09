#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@Author: Yang Bai
@Date: 2021.04.02
@Function: auto tests for all the input files in tests folder,
           all the release should be tested by this script!!!
"""
import os 
from pathlib import Path
import subprocess
import sys
import time

cpus=2
if len(sys.argv)>=3:
    if '-n' in sys.argv[2-1]:
        cpus=int(sys.argv[3-1])
        if cpus<1:
            cpus=1

currentdir=os.getcwd()
parrentdir=Path(currentdir).parent
if 'AsFem' not in str(parrentdir):
    parrentdir=currentdir
TestDir=str(parrentdir)+'/tests/'
AsFem=str(parrentdir)+'/bin/asfem'


print('**********************************************************************************')
print('*** We start to run the auto test script for all the test input file ...')
print('*** We are in folder:%s'%(currentdir))
print('*** Parent dir is: %s'%(parrentdir))
print('*** AsFem executable file is :%s'%(AsFem))
print('*** Test input file folder is :%s'%(TestDir))
print('*** Using %d cpus for auto-test'%(cpus))

timestart=time.time()

nFiles=0;nSucess=0
FailedFileList=[]
for subdir,dirs,files in os.walk(TestDir):
    print('***----------------------------------------------------------------------------***')
    print('***   start to run input files in %s'%(subdir))
    for file in files:
        if ('.i' in file) and ('.inp' not in file):
            arg=subdir+'/'+file
            print('***     running %s'%(file))
            os.chdir(subdir)
            if 'mesh' in subdir:
                args='mpirun -np %d '%(cpus)+AsFem+' -i '+file+' --read-only'
                result=subprocess.run(args,shell=True,capture_output=True) 
            else:
                args='mpirun -np %d '%(cpus)+AsFem+' -i '+file
                result=subprocess.run(args,shell=True,capture_output=True)
            nFiles+=1
            if ('AsFem exit due to some errors' in result.stdout.decode("utf-8")) or ('Error' in result.stdout.decode("utf-8")):
                sys.stdout.write("\033[1;31m") # set to red color
                print('***     %s is failed!'%(file))
                sys.stdout.write("\033[0;0m")  # reset color
                FailedFileList.append(file)
            else:
                nSucess+=1
                sys.stdout.write("\033[1;34m") # set to blue color
                print('***     %s is success!'%(file))
                sys.stdout.write("\033[0;0m")  # reset

timeend=time.time()
duration=timeend-timestart

print('**********************************************************************************')
print('*** Test finished, test files=%d, success=%d, failed=%d [elapse time=%13.5e]!'%(nFiles,nSucess,nFiles-nSucess,duration))
if len(FailedFileList)>0:
    print('*** The failed input files are:')
    print(FailedFileList)
print('**********************************************************************************')
