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

mpi="/home/by/Programs/openmpi/4.1.0/bin/mpirun"


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
                result=subprocess.run([AsFem,"-i",arg,"-read-only"],capture_output=True)
            else:
                result=subprocess.run([AsFem,"-i",arg],capture_output=True)
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

print('**********************************************************************************')
print('*** Test finished, test files=%d, sucess=%d, failed=%d !'%(nFiles,nSucess,nFiles-nSucess))
if len(FailedFileList)>0:
    print('*** The failed input files are:')
    print(FailedFileList)
print('**********************************************************************************')