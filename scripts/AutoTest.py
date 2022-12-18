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
scriptdir=Path(__file__).parent.resolve()
asfemdir=Path(scriptdir).parent
TestDir=str(asfemdir)+'/test_input/'
AsFem=str(asfemdir)+'/bin/asfem'


print('**********************************************************************************')
print('*** We start to run the auto test script for all the test input files ...')
print('*** We are in folder: %s'%(currentdir))
print('*** AutoTest script is in folder: %s'%(scriptdir))
print('*** AsFem dir is: %s'%(asfemdir))
print('*** AsFem executable file is : %s'%(AsFem))
print('*** Test input file folder is : %s'%(TestDir))
print('*** Using %d cpus for the auto-test'%(cpus))

timestart=time.time()

nFiles=0;nSucess=0
FailedFileList=[]
for subdir,dirs,files in os.walk(TestDir):
    print('***----------------------------------------------------------------------------***')
    print('***   start to run input files in %s'%(subdir))
    for file in files:
        if ('.json' in file) and ('input-template.json' not in file):
            arg=subdir+'/'+file
            print('***     running test for %s'%(file))
            os.chdir(subdir)
            if ('mesh' in subdir) or ('dofs' in subdir):
                args='mpirun -np %d '%(cpus)+AsFem+' -i '+file+' --read-only'
                result=subprocess.run(args,shell=True,capture_output=True) 
            else:
                args='mpirun -np %d '%(cpus)+AsFem+' -i '+file
                result=subprocess.run(args,shell=True,capture_output=True)
            nFiles+=1
            if ('AsFem exit due to some errors' in result.stdout.decode("utf-8")) or ('Error' in result.stdout.decode("utf-8")):
                sys.stdout.write("\033[1;31m") # set to red color
                print('***     %s fails !'%(file))
                sys.stdout.write("\033[0;0m")  # reset color
                FailedFileList.append(file)
            elif ('Static analysis is done' in result.stdout.decode("utf-8")) or ('Transient analysis is done' in result.stdout.decode("utf-8")) or ('\'Simulation is done\'' in result.stdout.decode("utf-8")) or ('AsFem has been executed in \'read-only\' mode' in result.stdout.decode("utf-8")):
                nSucess+=1
                sys.stdout.write("\033[1;34m") # set to blue color
                print('***     %s is done (success) !'%(file))
                sys.stdout.write("\033[0;0m")  # reset
            elif (len(result.stdout.decode("utf-8"))<2):
                # from cases with null message output
                nSucess+=1
                sys.stdout.write("\033[1;34m") # set to blue color
                print('***     %s is done (success) !'%(file))
                sys.stdout.write("\033[0;0m")  # reset
            else:
                sys.stdout.write("\033[1;31m") # set to red color
                print('***     %s fails (maybe no executable file?)!'%(file))
                print('*** The error message is:\n ',result.stdout.decode("utf-8"))
                sys.stdout.write("\033[0;0m")  # reset color
                FailedFileList.append(file)
                sys.exit()

timeend=time.time()
duration=timeend-timestart

print('**********************************************************************************')
print('*** Tests finished, test files=%d, success=%d, failure=%d [elapse time=%12.4e]!'%(nFiles,nSucess,nFiles-nSucess,duration))
if len(FailedFileList)>0:
    print('*** The failed input files are:')
    print(FailedFileList)
print('**********************************************************************************')
