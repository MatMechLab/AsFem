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
//+++ Date   : 2022.05.07
//+++ Purpose: the input file reading system in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

/**
 * for external packages and AsFem's headers
 */
#include "nlohmann/json.hpp"

#include "Utils/MessagePrinter.h"
#include "Utils/StringUtils.h"

#include "Mesh/Mesh.h"
#include "DofHandler/DofHandler.h"
#include "ElmtSystem/ElmtSystem.h"
#include "FE/FE.h"
#include "BCSystem/BCSystem.h"
#include "ICSystem/ICSystem.h"
#include "ProjectionSystem/ProjectionSystem.h"
#include "NonlinearSolver/NonlinearSolver.h"
#include "TimeStepping/TimeStepping.h"
#include "OutputSystem/OutputSystem.h"
#include "Postprocess/Postprocessor.h"
#include "FEProblem/FEJobBlock.h"

using std::string;
using std::vector;

/**
 * this class implement the input file reading for AsFem
 */
class InputSystem{
public:
    /**
     * constructor
     */
    InputSystem();
    InputSystem(int args,char *argv[]);
    /**
     * deconstructor
     */
    ~InputSystem();
    //*******************************************************
    //*** general settings
    //*******************************************************
    /**
     * initialize the input file reading system
     * @param args integer number of total argv
     * @param argv char vector taken from command line
     */
    void init(int args,char *argv[]);
    //*******************************************************
    //*** reading function for different blocks
    //*******************************************************
    /**
     * read the input file (json)
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofhandler class
     * @param t_elmtSystem the element system class
     * @param t_fe the fe space class for shape function and qpoints
     * @param t_bcsystem the boundary condition system class
     * @param t_icsystem the initial condition system class
     * @param t_projsystem the projection system class
     * @param t_nlsolver the nonlinear solver system
     * @param t_timestepping the time stepping system
     * @param t_output the output system
     * @param t_postprocess the postprocessor system
     * @param t_jobblock the job block
     */
    bool readInputFile(Mesh &t_mesh,DofHandler &t_dofhandler,ElmtSystem &t_elmtSystem,
                       FE &t_fe,
                       BCSystem &t_bcsystem,ICSystem &t_icsystem,
                       ProjectionSystem &t_projsystem,
                       NonlinearSolver &t_nlsolver,
                       TimeStepping &t_timestepping,
                       OutputSystem &t_output,
                       Postprocessor &t_postprocess,
                       FEJobBlock &t_jobblock);
    //*******************************************************
    //*** gettings
    //*******************************************************
    /**
     * get the name of input file, a string is returned
     */
    inline string getInputFileName()const{return m_inputfile_name;}

    /**
     * check whether the read-only flag is true
     */
    bool isReadOnly()const{return m_readonly;}

private:
    /**
     * read the mesh block from json file
     * @param t_json the json parse which contains 'mesh'
     * @param t_mesh the mesh class
     */
    bool readMeshBlock(nlohmann::json &t_json,Mesh &t_mesh);
    /**
     * read the dofs block from json file
     * @param t_json the json parse which contains 'dofs'
     * @param t_dofhandler the dofHandler class
     */
    bool readDofsBlock(nlohmann::json &t_json,DofHandler &t_dofhandler);
    /**
     * read the element block from json file
     * @param t_json the json parse which contains 'elements'
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofHandler class
     * @param t_elmtsystem the element system class
     */
    bool readElmtsBlock(nlohmann::json &t_json,const Mesh &t_mesh,const DofHandler &t_dofhandler,ElmtSystem &t_elmtsystem);
    /**
     * read the qpoint block
     * @param t_json the json parse which contains 'elements'
     * @param t_mesh the mesh class
     * @param t_fe the fe class
     */
    bool readQPointBlock(nlohmann::json &t_json,const Mesh &t_mesh,FE &t_fe);
    /**
     * read the shape function block
     * @param t_json the json parse which contains 'elements'
     * @param t_mesh the mesh class
     * @param t_fe the fe class
     */
    bool readShapeFunBlock(nlohmann::json &t_json,const Mesh &t_mesh,FE &t_fe);
    /**
     * read the boundary condition blocks
     * @param t_json the json parse which contains 'elements'
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofHandler class
     * @param t_bcsystem the boundary condition class
     */
    bool readBCsBlock(nlohmann::json &t_json,const Mesh &t_mesh,const DofHandler &t_dofhandler,BCSystem &t_bcsystem);
    /**
     * read the initial condition blocks
     * @param t_json the json parse which contains 'elements'
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofHandler class
     * @param t_icsystem the initial condition class
     */
    bool readICsBlock(nlohmann::json &t_json,const Mesh &t_mesh,const DofHandler &t_dofhandler,ICSystem &t_icsystem);
    /**
     * read the projection blocks
     * @param t_json the json parse which contains 'elements'
     * @param t_dofhandler the dofHandler class
     * @param t_projsystem the projection system class
     */
    bool readProjectionBlock(nlohmann::json &t_json,const DofHandler &t_dofhandler,ProjectionSystem &t_projsystem);
    /**
     * read the nonlinear solver blocks
     * @param t_json the json parse which contains 'elements'
     * @param t_nlsolver the nonlinear solver class
     */
    bool readNLSolverBlock(nlohmann::json &t_json,NonlinearSolver &t_nlsolver);
    /**
     * read the time stepping blocks
     * @param t_json the json parse which contains 'elements'
     * @param t_timestepping the time stepping solver class
     */
    bool readTimeSteppingBlock(nlohmann::json &t_json,TimeStepping &t_timestepping);
    /**
     * read the output solver blocks
     * @param t_json the json parse which contains 'elements'
     * @param t_output the output system class
     */
    bool readOutputBlock(nlohmann::json &t_json,OutputSystem &t_output);
    /**
     * read the postprocessor blocks
     * @param t_json the json parse which contains 'elements'
     * @param t_mesh the mesh class
     * @param t_dofhandler the dofhandler class
     * @param t_postprocessor the postprocessor class
     */
    bool readPostprocessBlock(nlohmann::json &t_json,const Mesh &t_mesh,const DofHandler &t_dofhandler,Postprocessor &t_postprocessor);
    /**
     * read the fe job block
     * @param t_json the json parse which contains 'elements'
     * @param t_jobblock the nonlinear solver class
     */
    bool readJobBlock(nlohmann::json &t_json,FEJobBlock &t_jobblock);

private:
    bool m_readonly;/**< boolean flag, if readonly=true, it will only read the mesh block */
    bool m_hasinputfile;/**< boolean flag, if true then the input file is loaded */
    string m_inputfile_name;/**< string for the name of input file */
    string m_meshfile_name;/**< string for the name of mesh file(external mesh file)*/
    nlohmann::json m_json;/**< json file reader */
    
};