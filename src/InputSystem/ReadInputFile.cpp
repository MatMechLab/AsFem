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
//+++ Purpose: read the input file
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readInputFile(Mesh &t_mesh,DofHandler &t_dofhandler,
                                ElmtSystem &t_elmtSystem,
                                FE &t_fe,
                                BCSystem &t_bcsystem,
                                ICSystem &t_icsystem,
                                ProjectionSystem &t_projsystem,
                                NonlinearSolver &t_nlsolver,
                                TimeStepping &t_timestepping,
                                OutputSystem &t_output,
                                Postprocessor &t_postprocess,
                                FEJobBlock &t_jobblock){
    ifstream in;
    if(m_hasinputfile){
        in.open(m_inputfile_name.c_str(),ios::in);
        while(!in.is_open()){
            MessagePrinter::printWarningTxt("can\'t open the input file(name="+m_inputfile_name+")");
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the correct input file name:");
            cin>>m_inputfile_name;
            in.open(m_inputfile_name.c_str(),ios::in);
        }
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the correct input file name:");
        cin>>m_inputfile_name;
        in.open(m_inputfile_name.c_str(),ios::in);
        while(!in.is_open()){
            MessagePrinter::printWarningTxt("can\'t open the input file");
            PetscPrintf(PETSC_COMM_WORLD,"*** Please enter the input file name:");
            cin>>m_inputfile_name;
        }
        m_hasinputfile=true;
    }

    string str;
    stringstream strStream;
    strStream<<in.rdbuf();
    str=strStream.str();
    in.close();
    m_json=nlohmann::json::parse(str);

    bool HasMeshBlock=false;
    bool HasDofsBlock=false;
    bool HasElmtsBlock=false;
    bool HasShapeFunBlock=false;
    bool HasQPointBlock=false;
    bool HasBCBlock=false;
    bool HasICBlock=false;
    bool HasProjectionBlock=false;
    bool HasNLSolverBlock=false;
    bool HasJobBlock=false;
    bool HasOutputBlock=false;
    bool HasPPSBlock=false;
    bool HasTimeSteppingBlock=false;

    // check whether all the blocks name are valid
    for(auto it=m_json.begin();it!=m_json.end();it++){
        if(it.key()=="mesh"||
           it.key()=="dofs"||
           it.key()=="elements"||
           it.key()=="qpoints"||
           it.key()=="shapefuns"||
           it.key()=="bcs"||
           it.key()=="ics"||
           it.key()=="projection"||
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


    if(m_json.contains("mesh")){
        // read the mesh block
        if(readMeshBlock(m_json.at("mesh"),t_mesh)){
            HasMeshBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'mesh' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("dofs")){
        // read the dofs block
        if(readDofsBlock(m_json.at("dofs"),t_dofhandler)){
            HasDofsBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'dofs' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("elements")){
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

        if(readElmtsBlock(m_json.at("elements"),t_mesh,t_dofhandler,t_elmtSystem)){
            HasElmtsBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'elements' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("qpoints")){
        // read the qpoints block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'qpoints' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(readQPointBlock(m_json.at("qpoints"),t_mesh,t_fe)){
            HasQPointBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'qpoints' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("shapefuns")){
        // read the shapefuns block
        if(!HasMeshBlock){
            MessagePrinter::printErrorTxt("you must define the 'mesh' block before your 'shapefuns' block, "
                                          "please check your input file");
            MessagePrinter::exitAsFem();
        }
        if(readShapeFunBlock(m_json.at("shapefuns"),t_mesh,t_fe)){
            HasShapeFunBlock=true;
        }
        else{
            MessagePrinter::printErrorTxt("something is incorrect in your 'shapefuns' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("bcs")){
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

        if(readBCsBlock(m_json.at("bcs"),t_mesh,t_dofhandler,t_bcsystem)){
            HasBCBlock=true;
        }
        else{
            HasBCBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'bcs' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("ics")){
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

        if(readICsBlock(m_json.at("ics"),t_mesh,t_dofhandler,t_icsystem)){
            HasBCBlock=true;
        }
        else{
            HasICBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'ics' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("projection")){
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

        if(readProjectionBlock(m_json.at("projection"),t_dofhandler,t_projsystem)){
            HasProjectionBlock=true;
        }
        else{
            HasProjectionBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'projection' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("nlsolver")){
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

        if(readNLSolverBlock(m_json.at("nlsolver"),t_nlsolver)){
            HasNLSolverBlock=true;
        }
        else{
            HasNLSolverBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'nlsolver' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }

    t_postprocess.setInputFileName(m_inputfile_name);
    if(m_json.contains("postprocess")){
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

        if(readPostprocessBlock(m_json.at("postprocess"),t_mesh,t_dofhandler,t_postprocess)){
            HasPPSBlock=true;
        }
        else{
            HasPPSBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'postprocess' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("timestepping")){
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

        if(readTimeSteppingBlock(m_json.at("timestepping"),t_timestepping)){
            HasTimeSteppingBlock=true;
        }
        else{
            HasTimeSteppingBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'timestepping' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("output")){
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

        if(readOutputBlock(m_json.at("output"),t_output)){
            HasOutputBlock=true;
        }
        else{
            HasOutputBlock=false;
            MessagePrinter::printErrorTxt("something is incorrect in your 'output' block, please check your input file");
            MessagePrinter::exitAsFem();
        }
    }
    
    if(m_json.contains("job")){
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

        if(readJobBlock(m_json.at("job"),t_jobblock)){
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
        t_fe.m_bulk_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        t_fe.m_bulk_qpoints.setDim(t_mesh.getBulkMeshMaxDim());
        t_fe.m_bulk_qpoints.setOrder(t_mesh.getBulkMeshBulkElmtOrder()+2);
        t_fe.m_bulk_qpoints.setMeshType(t_mesh.getBulkMeshBulkElmtMeshType());

        t_fe.m_surface_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        t_fe.m_surface_qpoints.setDim(2);
        t_fe.m_surface_qpoints.setOrder(t_mesh.getBulkMeshBulkElmtOrder()+2);
        t_fe.m_surface_qpoints.setMeshType(t_mesh.getBulkMeshSurfaceElmtMeshType());

        t_fe.m_line_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        t_fe.m_line_qpoints.setDim(1);
        t_fe.m_line_qpoints.setOrder(t_mesh.getBulkMeshBulkElmtOrder()+2);
        t_fe.m_line_qpoints.setMeshType(t_mesh.getBulkMeshLineElmtMeshType());
    }
    if(!HasShapeFunBlock){
        // use default value for shapefuns
        t_fe.m_bulk_shp.setMeshType(t_mesh.getBulkMeshBulkElmtMeshType());
        t_fe.m_bulk_shp.setShapeFunType(ShapeFunType::DEFAULT);

        t_fe.m_surface_shp.setMeshType(t_mesh.getBulkMeshSurfaceElmtMeshType());
        t_fe.m_surface_shp.setShapeFunType(ShapeFunType::DEFAULT);

        t_fe.m_line_shp.setMeshType(t_mesh.getBulkMeshLineElmtMeshType());
        t_fe.m_line_shp.setShapeFunType(ShapeFunType::DEFAULT);
    }

    if(!HasBCBlock){
        MessagePrinter::printWarningTxt("no [bcs] block found in your input file, then the 'zero' neumann bc is assumed");
    }

    if(!HasICBlock){
        MessagePrinter::printWarningTxt("no [ics] block found in your input file, then no any initial conditions will be applied");
    }

    if(!HasProjectionBlock){
        MessagePrinter::printWarningTxt("no [projection] block found in your input file, then no quantities will be projected");
    }

    if(!HasPPSBlock){
        MessagePrinter::printWarningTxt("no [postprocess] block found in your input file, then no postprocess will be executed");
    }

    if(!HasNLSolverBlock){
        MessagePrinter::printWarningTxt("no [nlsolver] block found in your input file, then the default options will be used");
        t_nlsolver.m_nlsolverblock.init();
    }

    if(!HasOutputBlock){
        MessagePrinter::printWarningTxt("no [output] block found in your input file, then the default options will be used");
        t_output.setFileFormat(ResultFileFormat::VTU);
        t_output.setIntervalNum(1);
    }
    t_output.setInputFileName(m_inputfile_name);

    if(t_jobblock.m_jobtype==FEJobType::TRANSIENT&&!HasTimeSteppingBlock){
        MessagePrinter::printErrorTxt("no [timestepping] block found in your input file, one can\'t do transient analysis without a timestepping block");
        MessagePrinter::exitAsFem();
    }

    if(!HasJobBlock && !m_readonly){
        MessagePrinter::printErrorTxt("no [job] block found in your input file, one can\'t do FEM analysis without a job");
        MessagePrinter::exitAsFem();
    }

    return HasMeshBlock;

}