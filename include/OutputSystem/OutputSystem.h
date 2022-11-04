//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group@CopyRight 2020-present
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2020.07.12
//+++ Purpose: define the output system for AsFem, where all the 
//+++          results should be written out to the result file
//+++          by this class
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include "OutputSystem/ResultFileFormat.h"
#include "OutputSystem/VTUWriter.h"


/**
 * This class implement the result output function
 */
class OutputSystem:public VTUWriter{
public:
    /**
     * constructor
     */
    OutputSystem();

    /**
     * setup the input file name
     * @param filename the string name of input file
     */
    void setInputFileName(const string &filename){m_inputfile_name=filename;}
    /**
     * setup the output interval
     */
    void setIntervalNum(const int &interval){m_intervals=interval;}
    /**
     * setup the output file format
     * @param t_format the file format type
     */
    void setFileFormat(const ResultFileFormat &t_format){m_fileformat=t_format;}

    /**
     * save result to different files according to the output format
     * @param t_step the time step, -1 for static analysis output
     * @param t_mesh the mesh class
     * @param t_dofHandler the dofhandler class
     * @param t_solution the solution class
     * @param t_projection the projection class
     */
    void saveResults2File(const int &t_step,
                          const Mesh &t_mesh,
                          const DofHandler &t_dofHandler,
                          SolutionSystem &t_solution,
                          ProjectionSystem &t_projection);
    
    /**
     * save result info to pvd file, this is used by vtu file
     * @param current_time the time of current time step
     */
    void savePVDResults(const double &current_time);
    /**
     * save the header info of pvd file
     */
    void savePVDHead();
    /**
     * save the end info of pvd file
     */
    void savePVDEnd();

    //*****************************************************
    //*** general gettings
    //*****************************************************
    /**
     * get the string name of the output file
     */
    inline string getOutputFileName()const{return m_outputfile_name;}
    /**
     * get the output interval number
     */
    inline int getIntervalNum()const{return m_intervals;}

    /**
     * print out the output system info
     */
    void printInfo()const;


private:
    string m_inputfile_name;/**< for the string name of input file */
    string m_outputfile_name;/**< for the string name of output file */
    string m_pvdfile_name;/**< for the stringname of pvd file */
    ResultFileFormat m_fileformat;/**< for the result file format */

    int m_intervals;/**< the output interval number */

private:
    PetscMPIInt m_rank;/** for the processor id */

};