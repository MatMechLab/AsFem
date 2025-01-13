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
//+++ Date   : 2022.05.07
//+++ Purpose: read the input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"
#include "MPIUtils/MPIDataBus.h"

bool InputSystem::readInputFile(FECell &t_fecell,DofHandler &t_dofhandler,
                                ElmtSystem &t_elmtSystem,
                                FE &t_fe,
                                BCSystem &t_bcsystem,
                                ICSystem &t_icsystem,
                                ProjectionSystem &t_projsystem,
                                LinearSolver &t_linearsolver,
                                NonlinearSolver &t_nlsolver,
                                TimeStepping &t_timestepping,
                                OutputSystem &t_output,
                                Postprocessor &t_postprocess,
                                FEJobBlock &t_jobblock){
    ifstream in;
    if(m_HasInputFile){
        in.open(m_InputFileName.c_str(),ios::in);
        while(!in.is_open()){
            MessagePrinter::printWarningTxt("can\'t open the input file(name="+m_InputFileName+")");
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the correct input file name:");
            cin>>m_InputFileName;
            in.open(m_InputFileName.c_str(),ios::in);
        }
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the correct input file name:");
        cin>>m_InputFileName;
        in.open(m_InputFileName.c_str(),ios::in);
        while(!in.is_open()){
            MessagePrinter::printWarningTxt("can\'t open the input file");
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the input file name:");
            cin>>m_InputFileName;
        }
        m_HasInputFile=true;
    }


    string str;
    stringstream strStream;
    strStream<<in.rdbuf();
    str=strStream.str();
    in.close();
    m_Json=nlohmann::json::parse(str);

    bool HasMeshBlock=false;
    bool HasDofsBlock=false;
    bool HasElmtsBlock=false;
    bool HasShapeFunBlock=false;
    bool HasQPointBlock=false;
    bool HasBCBlock=false;
    bool HasICBlock=false;
    bool HasProjectionBlock=false;
    bool HasLinearSolverBlock=false;
    bool HasNLSolverBlock=false;
    bool HasJobBlock=false;
    bool HasOutputBlock=false;
    bool HasPPSBlock=false;
    bool HasTimeSteppingBlock=false;

    // check whether all the blocks name are valid
    for(auto it=m_Json.begin();it!=m_Json.end();it++){
        if(it.key()=="mesh"||
           it.key()=="dofs"||
           it.key()=="elements"||
           it.key()=="qpoints"||
           it.key()=="shapefuns"||
           it.key()=="bcs"||
           it.key()=="ics"||
           it.key()=="projection"||
           it.key()=="linearsolver"||
           it.key()=="nlsolver"||
           it.key()=="timestepping"||
           it.key()=="output"||
           it.key()=="postprocess"||
           it.key()=="job"||
           it.key()=="description"){
            continue;
        }
        else{
            MessagePrinter::printErrorTxt("'"+it.key()+"' is not a valid block name, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }


    if(m_Json.contains("mesh")){
        // read the mesh block
        if(readMeshBlock(m_Json.at("mesh"),t_fecell)){
            HasMeshBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'mesh' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("dofs")){
        // read the dofs block
        if(readDofsBlock(m_Json.at("dofs"),t_dofhandler)){
            HasDofsBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'dofs' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("elements")){
        // read the elements block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'elements' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'elements' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(readElmtsBlock(m_Json.at("elements"),t_fecell,t_dofhandler,t_elmtSystem)){
            HasElmtsBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'elements' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("qpoints")){
        // read the qpoints block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'qpoints' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(readQPointBlock(m_Json.at("qpoints"),t_fecell,t_fe)){
            HasQPointBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'qpoints' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("shapefuns")){
        // read the shapefuns block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'shapefuns' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(readShapeFunBlock(m_Json.at("shapefuns"),t_fecell,t_fe)){
            HasShapeFunBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'shapefuns' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("bcs")){
        // read the bcs block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'bcs' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'bcs' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(readBCsBlock(m_Json.at("bcs"),t_fecell,t_dofhandler,t_bcsystem)){
            HasBCBlock=true;
        }
        else{
            HasBCBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'bcs' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("ics")){
        // read the ics block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'ics' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'ics' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(readICsBlock(m_Json.at("ics"),t_fecell,t_dofhandler,t_icsystem)){
            HasICBlock=true;
        }
        else{
            HasICBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'ics' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("projection")){
        // read the bcs block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'projection' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'projection' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(readProjectionBlock(m_Json.at("projection"),t_dofhandler,t_projsystem)){
            HasProjectionBlock=true;
        }
        else{
            HasProjectionBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'projection' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }

    if(m_Json.contains("linearsolver")){
        // read the bcs block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'linearsolver' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'linearsolver' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(readLinearSolverBlock(m_Json.at("linearsolver"),t_linearsolver)){
            HasLinearSolverBlock=true;
        }
        else{
            HasLinearSolverBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'linearsolver' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("nlsolver")){
        // read the bcs block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'nlsolver' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'nlsolver' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(readNLSolverBlock(m_Json.at("nlsolver"),t_nlsolver)){
            HasNLSolverBlock=true;
        }
        else{
            HasNLSolverBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'nlsolver' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }

    t_postprocess.setInputFileName(m_InputFileName);
    if(m_Json.contains("postprocess")){
        // read the bcs block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'postprocess' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'postprocess' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(readPostprocessBlock(m_Json.at("postprocess"),t_fecell,t_dofhandler,t_postprocess)){
            HasPPSBlock=true;
        }
        else{
            HasPPSBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'postprocess' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("timestepping")){
        // read the bcs block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'timestepping' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'timestepping' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(readTimeSteppingBlock(m_Json.at("timestepping"),t_timestepping)){
            HasTimeSteppingBlock=true;
        }
        else{
            HasTimeSteppingBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'timestepping' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("output")){
        // read the bcs block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'output' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'output' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(!HasElmtsBlock){
            MessagePrinter::printErrorTxt("you must define the 'elements' block before your 'output' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(readOutputBlock(m_Json.at("output"),t_output)){
            HasOutputBlock=true;
        }
        else{
            HasOutputBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'output' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_Json.contains("job")){
        // read the bcs block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'job' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(!HasDofsBlock){
            MessagePrinter::printErrorTxt("you must define the 'dofs' block before your 'job' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(!HasElmtsBlock){
            MessagePrinter::printErrorTxt("you must define the 'elements' block before your 'job' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(readJobBlock(m_Json.at("job"),t_jobblock)){
            HasJobBlock=true;
        }
        else{
            HasJobBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'job' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(HasDofsBlock){}
    if(HasElmtsBlock){}

    if(!HasQPointBlock){
        // use default options for qpoints
        t_fe.m_BulkQpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        t_fe.m_BulkQpoints.setDim(t_fecell.getFECellMaxDim());
        t_fe.m_BulkQpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+2);
        t_fe.m_BulkQpoints.setMeshType(t_fecell.getFECellBulkElmtMeshType());

        t_fe.m_SurfaceQpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        t_fe.m_SurfaceQpoints.setDim(2);
        t_fe.m_SurfaceQpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+2);
        t_fe.m_SurfaceQpoints.setMeshType(t_fecell.getFECellSurfElmtMeshType());

        t_fe.m_LineQpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        t_fe.m_LineQpoints.setDim(1);
        t_fe.m_LineQpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+2);
        t_fe.m_LineQpoints.setMeshType(t_fecell.getFECellLineElmtMeshType());
    }
    if(!HasShapeFunBlock){
        // use default value for shapefuns
        t_fe.m_BulkShp.setMeshType(t_fecell.getFECellBulkElmtMeshType());
        t_fe.m_BulkShp.setShapeFunType(ShapeFunType::DEFAULT);

        t_fe.m_SurfaceShp.setMeshType(t_fecell.getFECellSurfElmtMeshType());
        t_fe.m_SurfaceShp.setShapeFunType(ShapeFunType::DEFAULT);

        t_fe.m_LineShp.setMeshType(t_fecell.getFECellLineElmtMeshType());
        t_fe.m_LineShp.setShapeFunType(ShapeFunType::DEFAULT);
    }

    if(!HasBCBlock){
        MessagePrinter::printWarningTxt("no 'bcs' block found in your input file, then the 'zero' neumann bc will be applied");
    }

    if(!HasICBlock){
        MessagePrinter::printWarningTxt("no 'ics' block found in your input file, then no any initial conditions will be applied");
    }

    if(!HasProjectionBlock){
        MessagePrinter::printWarningTxt("no 'projection' block found in your input file, then no quantities will be projected");
    }

    if(!HasPPSBlock){
        MessagePrinter::printWarningTxt("no 'postprocess' block found in your input file, then no postprocess will be executed");
    }


    if (!HasLinearSolverBlock) {
        MessagePrinter::printWarningTxt("no 'linearsolver' block found in your input file, then the default options will be used");
        t_linearsolver.setDefaultParams();
    }

    if(!HasNLSolverBlock){
        MessagePrinter::printWarningTxt("no 'nlsolver' block found in your input file, then the default options will be used");
        t_nlsolver.m_NlSolverBlock.init();
    }

    if(!HasOutputBlock){
        MessagePrinter::printWarningTxt("no 'output' block found in your input file, then the default options will be used");
        t_output.setFileFormat(ResultFileFormat::VTU);
        t_output.setIntervalNum(1);
    }
    t_output.setInputFileName(m_InputFileName);

    if(t_jobblock.m_JobType==FEJobType::TRANSIENT&&!HasTimeSteppingBlock){
        MessagePrinter::printErrorTxt("no 'timestepping' block found in your input file, one can\'t do transient analysis without a timestepping block");
        MessagePrinter::exitAsFem();
    }

    if(t_jobblock.m_JobType==FEJobType::STATIC&&HasTimeSteppingBlock){
        MessagePrinter::printErrorTxt("'timestepping' block is found in your input file for a static job, this dosen\'t make sense");
        MessagePrinter::exitAsFem();
    }

    if(!HasJobBlock && !m_ReadOnly){
        MessagePrinter::printErrorTxt("no 'job' block found in your input file, one can\'t do FEM analysis without a job");
        MessagePrinter::exitAsFem();
    }

    return HasMeshBlock;

}