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
//+++ Date   : 2022.06.13
//+++ Purpose: read the qpoint block
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"

bool InputSystem::readQPointBlock(nlohmann::json &t_json,const FECell &t_fecell,FE &t_fe){
    // here json already contains "qpoints"

    // once user defines the "qpoints" block, we init all the sub-block with default options
    t_fe.m_bulk_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
    t_fe.m_bulk_qpoints.setDim(t_fecell.getFECellMaxDim());
    t_fe.m_bulk_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+2);
    t_fe.m_bulk_qpoints.setMeshType(t_fecell.getFECellBulkElmtMeshType());

    t_fe.m_surface_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
    t_fe.m_surface_qpoints.setDim(2);
    t_fe.m_surface_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+2);
    t_fe.m_surface_qpoints.setMeshType(t_fecell.getFECellSurfElmtMeshType());

    t_fe.m_line_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
    t_fe.m_line_qpoints.setDim(1);
    t_fe.m_line_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+2);
    t_fe.m_line_qpoints.setMeshType(t_fecell.getFECellLineElmtMeshType());

    if(t_json.contains("bulk")){
        auto json_bulk=t_json.at("bulk");
        if(!json_bulk.contains("type")){
            MessagePrinter::printErrorTxt("can\'t find type in your 'bulk' qpoint block, please check your input file");
            MessagePrinter::exitAsFem();
        }

        string qtype=json_bulk.at("type");
        if(qtype.find("gauss-legendre")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        }
        else if(qtype.find("gauss-lobatto")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::GAUSSLOBATTO);
        }
        else if(qtype.find("user1")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER1);
        }
        else if(qtype.find("user2")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER2);
        }
        else if(qtype.find("user3")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER3);
        }
        else if(qtype.find("user4")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER4);
        }
        else if(qtype.find("user5")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER5);
        }
        else{
            MessagePrinter::printErrorTxt("unsupported gauss point type(="+qtype+") for bulk mesh, please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(json_bulk.contains("order")){
            if(!json_bulk.at("order").is_number_integer()){
                MessagePrinter::printErrorTxt("the value of 'order=' option is invalid for bulk mesh, please check your input file");
                MessagePrinter::exitAsFem();
            }
            int order=json_bulk.at("order");
            if(order<0){
                MessagePrinter::printErrorTxt("order="+to_string(order)+" is invalid for gauss point order of bulk mesh, please check your input file");
                MessagePrinter::exitAsFem();
            }
            t_fe.m_bulk_qpoints.setOrder(order);
        }
        else{
            t_fe.m_bulk_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+1);
        }
    }

    if(t_json.contains("surface")){
        if(t_fecell.getFECellMaxDim()<3){
            MessagePrinter::printErrorTxt("the max dim of your mesh is "+to_string(t_fecell.getFECellMaxDim())+", you can\'t define qpoint for a surface mesh,"
                                          " please check your input file");
            MessagePrinter::exitAsFem();
        }

        auto json_surface=t_json.at("surface");
        if(!json_surface.contains("type")){
            MessagePrinter::printErrorTxt("can\'t find type in your 'surface' qpoint block, please check your input file");
            MessagePrinter::exitAsFem();
        }

        string qtype=json_surface.at("type");
        if(qtype.find("gauss-legendre")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        }
        else if(qtype.find("gauss-lobatto")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::GAUSSLOBATTO);
        }
        else if(qtype.find("user1")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER1);
        }
        else if(qtype.find("user2")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER2);
        }
        else if(qtype.find("user3")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER3);
        }
        else if(qtype.find("user4")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER4);
        }
        else if(qtype.find("user5")!=string::npos){
            t_fe.m_bulk_qpoints.setQPointType(QPointType::USER5);
        }
        else{
            MessagePrinter::printErrorTxt("unsupported gauss point type(="+qtype+") for surface mesh, please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(json_surface.contains("order")){
            if(!json_surface.at("order").is_number_integer()){
                MessagePrinter::printErrorTxt("the value of 'order=' option is invalid for surface mesh, please check your input file");
                MessagePrinter::exitAsFem();
            }
            int order=json_surface.at("order");
            if(order<0){
                MessagePrinter::printErrorTxt("order="+to_string(order)+" is invalid for gauss point order of surface mesh, please check your input file");
                MessagePrinter::exitAsFem();
            }
            t_fe.m_surface_qpoints.setOrder(order);
        }
        else{
            t_fe.m_surface_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+1);
        }

    }

    if(t_json.contains("line")){
        if(t_fecell.getFECellMaxDim()<2){
            MessagePrinter::printErrorTxt("the max dim of your mesh is "+to_string(t_fecell.getFECellMaxDim())+", you can\'t define qpoint for a line mesh,"
                                          " please check your input file");
            MessagePrinter::exitAsFem();
        }

        auto json_line=t_json.at("line");
        if(!json_line.contains("type")){
            MessagePrinter::printErrorTxt("can\'t find type in your 'line' qpoint block, please check your input file");
            MessagePrinter::exitAsFem();
        }

        string qtype=json_line.at("type");
        if(qtype.find("gauss-legendre")!=string::npos){
            t_fe.m_line_qpoints.setQPointType(QPointType::GAUSSLEGENDRE);
        }
        else if(qtype.find("gauss-lobatto")!=string::npos){
            t_fe.m_line_qpoints.setQPointType(QPointType::GAUSSLOBATTO);
        }
        else if(qtype.find("user1")!=string::npos){
            t_fe.m_line_qpoints.setQPointType(QPointType::USER1);
        }
        else if(qtype.find("user2")!=string::npos){
            t_fe.m_line_qpoints.setQPointType(QPointType::USER2);
        }
        else if(qtype.find("user3")!=string::npos){
            t_fe.m_line_qpoints.setQPointType(QPointType::USER3);
        }
        else if(qtype.find("user4")!=string::npos){
            t_fe.m_line_qpoints.setQPointType(QPointType::USER4);
        }
        else if(qtype.find("user5")!=string::npos){
            t_fe.m_line_qpoints.setQPointType(QPointType::USER5);
        }
        else{
            MessagePrinter::printErrorTxt("unsupported gauss point type(="+qtype+") for line mesh, please check your input file");
            MessagePrinter::exitAsFem();
        }

        if(json_line.contains("order")){
            if(!json_line.at("order").is_number_integer()){
                MessagePrinter::printErrorTxt("the value of 'order=' option is invalid for line mesh, please check your input file");
                MessagePrinter::exitAsFem();
            }
            int order=json_line.at("order");
            if(order<0){
                MessagePrinter::printErrorTxt("order="+to_string(order)+" is invalid for gauss point order of line mesh, please check your input file");
                MessagePrinter::exitAsFem();
            }
            t_fe.m_line_qpoints.setOrder(order);
        }
        else{
            t_fe.m_line_qpoints.setOrder(t_fecell.getFECellBulkMeshOrder()+1);
        }

    }
    return true;
}