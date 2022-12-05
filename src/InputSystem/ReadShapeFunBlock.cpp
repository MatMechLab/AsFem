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
//+++ Date   : 2022.06.13
//+++ Purpose: read the shapefunction block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readShapeFunBlock(nlohmann::json &t_json,const Mesh &t_mesh,FE &t_fe){
    // json already constains "shapefuns"
    // before we read the info, we should preset the default value for each "shp"
    t_fe.m_bulk_shp.setMeshType(t_mesh.getBulkMeshBulkElmtMeshType());
    t_fe.m_bulk_shp.setShapeFunType(ShapeFunType::DEFAULT);

    t_fe.m_surface_shp.setMeshType(t_mesh.getBulkMeshSurfaceElmtMeshType());
    t_fe.m_surface_shp.setShapeFunType(ShapeFunType::DEFAULT);

    t_fe.m_line_shp.setMeshType(t_mesh.getBulkMeshLineElmtMeshType());
    t_fe.m_line_shp.setShapeFunType(ShapeFunType::DEFAULT);

    if(t_json.contains("bulk")){
        auto json_bulk=t_json.at("bulk");
        if(!json_bulk.contains("type")){
            MessagePrinter::printErrorTxt("can\'t find 'type=' option for your 'bulk' shapefun, please check your input file");
            MessagePrinter::exitAsFem();
        }
        string shptype=json_bulk.at("type");
        if(shptype.find("default")!=string::npos){
            t_fe.m_bulk_shp.setShapeFunType(ShapeFunType::DEFAULT);
        }
        else if(shptype.find("user1")!=string::npos){
            t_fe.m_bulk_shp.setShapeFunType(ShapeFunType::USER1);
        }
        else if(shptype.find("user2")!=string::npos){
            t_fe.m_bulk_shp.setShapeFunType(ShapeFunType::USER2);
        }
        else if(shptype.find("user3")!=string::npos){
            t_fe.m_bulk_shp.setShapeFunType(ShapeFunType::USER3);
        }
        else if(shptype.find("user4")!=string::npos){
            t_fe.m_bulk_shp.setShapeFunType(ShapeFunType::USER4);
        }
        else if(shptype.find("user5")!=string::npos){
            t_fe.m_bulk_shp.setShapeFunType(ShapeFunType::USER5);
        }
        else{
            MessagePrinter::printErrorTxt("unsupported shape function type(="+shptype+") in your 'bulk' shapefun, please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(json_bulk.contains("functions")){
            if(shptype.find("default")!=string::npos){
                MessagePrinter::printErrorTxt("for default shapefunction, you don\'t need to define the functions number"
                                              ", please check your input file");
                MessagePrinter::exitAsFem();
            }
            else{
                if(!json_bulk.at("functions").is_number_integer()){
                    MessagePrinter::printErrorTxt("invalid integer for functions in your 'bulk' shapefun, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                int funs=json_bulk.at("functions");
                if(funs<1){
                    MessagePrinter::printErrorTxt("functions="+to_string(funs)+" is invalid for your 'bulk' shapefun, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                t_fe.m_bulk_shp.setShapeFunNums(funs);
            }
        }
        else{
            if(shptype.find("user1")!=string::npos||
               shptype.find("user2")!=string::npos||
               shptype.find("user3")!=string::npos||
               shptype.find("user4")!=string::npos||
               shptype.find("user5")!=string::npos){
                MessagePrinter::printErrorTxt("can\'t find 'functions' in your 'bulk' shapefun, "
                                              "the user-defined shapefun must have 'functions'," 
                                              " please check your input file");
                MessagePrinter::exitAsFem();
            }
        }
    }

    if(t_json.contains("surface")){
        if(t_mesh.getBulkMeshMaxDim()<3){
            MessagePrinter::printErrorTxt("you can\'t define surface shape fun, the max dim of your mesh is "
                                          +to_string(t_mesh.getBulkMeshMaxDim())
                                          +", please check your input file");
            MessagePrinter::exitAsFem();
        }
        auto json_surface=t_json.at("surface");
        if(!json_surface.contains("type")){
            MessagePrinter::printErrorTxt("can\'t find 'type=' option for your 'surface' shapefun, please check your input file");
            MessagePrinter::exitAsFem();
        }
        string shptype=json_surface.at("type");
        if(shptype.find("default")!=string::npos){
            t_fe.m_surface_shp.setShapeFunType(ShapeFunType::DEFAULT);
        }
        else if(shptype.find("user1")!=string::npos){
            t_fe.m_surface_shp.setShapeFunType(ShapeFunType::USER1);
        }
        else if(shptype.find("user2")!=string::npos){
            t_fe.m_surface_shp.setShapeFunType(ShapeFunType::USER2);
        }
        else if(shptype.find("user3")!=string::npos){
            t_fe.m_surface_shp.setShapeFunType(ShapeFunType::USER3);
        }
        else if(shptype.find("user4")!=string::npos){
            t_fe.m_surface_shp.setShapeFunType(ShapeFunType::USER4);
        }
        else if(shptype.find("user5")!=string::npos){
            t_fe.m_surface_shp.setShapeFunType(ShapeFunType::USER5);
        }
        else{
            MessagePrinter::printErrorTxt("unsupported shape function type(="+shptype+") in your 'surface' shapefun, please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(json_surface.contains("functions")){
            if(shptype.find("default")!=string::npos){
                MessagePrinter::printErrorTxt("for default shapefunction, you don\'t need to define the functions number"
                                              ", please check your input file");
                MessagePrinter::exitAsFem();
            }
            else{
                if(!json_surface.at("functions").is_number_integer()){
                    MessagePrinter::printErrorTxt("invalid integer for functions in your 'surface' shapefun, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                int funs=json_surface.at("functions");
                if(funs<1){
                    MessagePrinter::printErrorTxt("functions="+to_string(funs)+" is invalid for your 'surface' shapefun, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                t_fe.m_surface_shp.setShapeFunNums(funs);
            }
        }
        else{
            if(shptype.find("user1")!=string::npos||
               shptype.find("user2")!=string::npos||
               shptype.find("user3")!=string::npos||
               shptype.find("user4")!=string::npos||
               shptype.find("user5")!=string::npos){
                MessagePrinter::printErrorTxt("can\'t find 'functions' in your 'surface' shapefun, "
                                              "the user-defined shapefun must have 'functions'," 
                                              " please check your input file");
                MessagePrinter::exitAsFem();
            }
        }
    }

    if(t_json.contains("line")){
        if(t_mesh.getBulkMeshMaxDim()<2){
            MessagePrinter::printErrorTxt("you can\'t define line shape fun, the max dim of your mesh is "
                                          +to_string(t_mesh.getBulkMeshMaxDim())
                                          +", please check your input file");
            MessagePrinter::exitAsFem();
        }
        auto json_line=t_json.at("line");
        if(!json_line.contains("type")){
            MessagePrinter::printErrorTxt("can\'t find 'type=' option for your 'line' shapefun, please check your input file");
            MessagePrinter::exitAsFem();
        }
        string shptype=json_line.at("type");
        if(shptype.find("default")!=string::npos){
            t_fe.m_line_shp.setShapeFunType(ShapeFunType::DEFAULT);
        }
        else if(shptype.find("user1")!=string::npos){
            t_fe.m_line_shp.setShapeFunType(ShapeFunType::USER1);
        }
        else if(shptype.find("user2")!=string::npos){
            t_fe.m_line_shp.setShapeFunType(ShapeFunType::USER2);
        }
        else if(shptype.find("user3")!=string::npos){
            t_fe.m_line_shp.setShapeFunType(ShapeFunType::USER3);
        }
        else if(shptype.find("user4")!=string::npos){
            t_fe.m_line_shp.setShapeFunType(ShapeFunType::USER4);
        }
        else if(shptype.find("user5")!=string::npos){
            t_fe.m_line_shp.setShapeFunType(ShapeFunType::USER5);
        }
        else{
            MessagePrinter::printErrorTxt("unsupported shape function type(="+shptype+") in your 'line' shapefun, please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(json_line.contains("functions")){
            if(shptype.find("default")!=string::npos){
                MessagePrinter::printErrorTxt("for default shapefunction, you don\'t need to define the functions number"
                                              ", please check your input file");
                MessagePrinter::exitAsFem();
            }
            else{
                if(!json_line.at("functions").is_number_integer()){
                    MessagePrinter::printErrorTxt("invalid integer for functions in your 'line' shapefun, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                int funs=json_line.at("functions");
                if(funs<1){
                    MessagePrinter::printErrorTxt("functions="+to_string(funs)+" is invalid for your 'line' shapefun, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                t_fe.m_line_shp.setShapeFunNums(funs);
            }
        }
        else{
            if(shptype.find("user1")!=string::npos||
               shptype.find("user2")!=string::npos||
               shptype.find("user3")!=string::npos||
               shptype.find("user4")!=string::npos||
               shptype.find("user5")!=string::npos){
                MessagePrinter::printErrorTxt("can\'t find 'functions' in your 'line' shapefun, "
                                              "the user-defined shapefun must have 'functions'," 
                                              " please check your input file");
                MessagePrinter::exitAsFem();
            }
        }
    }

    return true;
}