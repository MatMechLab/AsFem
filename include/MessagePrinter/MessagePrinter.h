//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai @ CopyRight 2020
//* https://github.com/yangbai90/AsFem.git
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************

#ifndef ASFEM_MESSAGEPRINTER_H
#define ASFEM_MESSAGEPRINTER_H

#include <iostream>
#include "petsc.h"

using namespace std;

void Msg_AsFem_Exit();


//*******************************************
//*** For input file
//*******************************************
void Msg_Input_InvalidArgs();
void Msg_Input_InvalidInputFileName();
void Msg_InputFileName_Invalid(string filename);

void Msg_Input_LineError(const PetscInt &linenum);

void Msg_Input_BlockBracketNotComplete();


#endif // ASFEM_MESSAGEPRINTER_H