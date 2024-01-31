//****************************************************************
//* This file is part of the AsFem framework
//* Advanced Simulation kit based on Finite Element Method (AsFem)
//* All rights reserved, Yang Bai/MM-Lab@CopyRight 2020-present
//* https://github.com/MatMechLab/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2024.01.31
//+++ Purpose: the application level management for AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once


#include "petsc.h"
#include "Utils/MessagePrinter.h"

/**
 * This class manipulate the app level process, i.e. init and finalize the program
*/
class Application{
public:
    /**
     * Constructor
    */
    Application();

    /**
     * Initialize the asfem app
     * @param args number of arguments
     * @param argv the argument string vector
    */
    PetscErrorCode init(int args,char *argv[]);
    /**
     * Finalize the asfem app
    */
    PetscErrorCode finalize();


    /**
     * print the welcome message and app info in your terminal
     * @param year the integer number of release year
     * @param month the integer number of release month
     * @param day the integer number of release day
     * @param version the double number of release version (2-digital)
    */
    void printAppInfo(const int &year,const int &month,const int &day,const double &version)const;

private:
    PetscOptions m_Options;/**< the petsc options structure */
};