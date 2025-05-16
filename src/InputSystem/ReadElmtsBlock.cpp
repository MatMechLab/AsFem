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
//+++ Date   : 2022.05.12
//+++ Purpose: read the elements block from json
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"
#include "ElmtSystem/ElmtBlock.h"

bool InputSystem::readElmtsBlock(nlohmann::json &t_json,const FECell &t_fecell,const DofHandler &t_dofhandler,ElmtSystem &t_elmtsystem){
    // now the json should already read 'elements'
    ElmtBlock elmtBlock;
    int blocks=0;
    bool HasType=false;
    bool HasDofs=false;
    for(auto it=t_json.begin();it!=t_json.end();it++){
        elmtBlock.reset();
        blocks+=1;

        HasType=false;
        HasDofs=false;
        
        elmtBlock.m_ElmtBlockIndex=blocks;
        elmtBlock.m_ElmtBlockName=it.key();
        nlohmann::json ejson=t_json.at(elmtBlock.m_ElmtBlockName);//it.value();// now ejson gose into 'elmt-x'
        if(ejson.contains("type")){
            if(!ejson.at("type").is_string()){
                MessagePrinter::printErrorTxt("invalid type name in element block-"+to_string(blocks)+", please check your input file");
                MessagePrinter::exitAsFem();
            }
            elmtBlock.m_ElmtTypeName=ejson.at("type");
            elmtBlock.m_ElmtType=ElmtType::NULLELMT;
            if(elmtBlock.m_ElmtTypeName=="poisson"){
                elmtBlock.m_ElmtType=ElmtType::POISSONELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="diffusion"){
                elmtBlock.m_ElmtType=ElmtType::DIFFUSIONELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="thermal"){
                elmtBlock.m_ElmtType=ElmtType::THERMALELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="mechanics"){
                elmtBlock.m_ElmtType=ElmtType::MECHANICSELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="dynamics"){
                elmtBlock.m_ElmtType=ElmtType::DYNAMICMECHANICSELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="allencahn"){
                elmtBlock.m_ElmtType=ElmtType::ALLENCAHNELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="cahnhilliard"){
                elmtBlock.m_ElmtType=ElmtType::CAHNHILLIARDELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="wave"){
                elmtBlock.m_ElmtType=ElmtType::WAVEELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="kobayashi"){
                elmtBlock.m_ElmtType=ElmtType::KOBAYASHIELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="laplace"){
                elmtBlock.m_ElmtType=ElmtType::LAPLACEELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="scalarbodysource"){
                elmtBlock.m_ElmtType=ElmtType::SCALARBODYSOURCEELMT;
            }
            // for coupled cases
            else if(elmtBlock.m_ElmtTypeName=="stresscahnhilliard"){
                elmtBlock.m_ElmtType=ElmtType::STRESSCAHNHILLIARDELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="stressdiffusion"){
                elmtBlock.m_ElmtType=ElmtType::STRESSDIFFUSIONELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="allencahnfracture"){
                elmtBlock.m_ElmtType=ElmtType::ALLENCAHNFRACTUREELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="miehefracture"){
                elmtBlock.m_ElmtType=ElmtType::MIEHEFRACTUREELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="diffusionacfracture"){
                elmtBlock.m_ElmtType=ElmtType::DIFFUSIONACFRACTUREELMT;
            }
            // for user-defined-element
            else if(elmtBlock.m_ElmtTypeName=="user1"){
                elmtBlock.m_ElmtType=ElmtType::USER1ELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="user2"){
                elmtBlock.m_ElmtType=ElmtType::USER2ELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="user3"){
                elmtBlock.m_ElmtType=ElmtType::USER3ELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="user4"){
                elmtBlock.m_ElmtType=ElmtType::USER4ELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="user5"){
                elmtBlock.m_ElmtType=ElmtType::USER5ELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="user6"){
                elmtBlock.m_ElmtType=ElmtType::USER6ELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="user7"){
                elmtBlock.m_ElmtType=ElmtType::USER7ELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="user8"){
                elmtBlock.m_ElmtType=ElmtType::USER8ELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="user9"){
                elmtBlock.m_ElmtType=ElmtType::USER9ELMT;
            }
            else if(elmtBlock.m_ElmtTypeName=="user10"){
                elmtBlock.m_ElmtType=ElmtType::USER10ELMT;
            }
            else{
                MessagePrinter::printErrorTxt("unsupported element type in "+elmtBlock.m_ElmtBlockName+", please check your input file");
                MessagePrinter::exitAsFem();
            }//end-of-element-type-name
            HasType=true;
        }
        else{
            MessagePrinter::printErrorTxt("can\'t find type name in element block-"+to_string(blocks)+", please check your input file");
            MessagePrinter::exitAsFem();
        }//end-of-'type'-reading

        if(ejson.contains("dofs")){
            if(ejson.at("dofs").size()<1){
                MessagePrinter::printErrorTxt("invalid dofs name or empty dofs name in element block-"+to_string(blocks)+", please check your input file");
                MessagePrinter::exitAsFem();
            }
            else{
                string dofname;
                for(int i=0;i<static_cast<int>(ejson.at("dofs").size());i++){
                    if(!ejson.at("dofs").at(i).is_string()){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of elmt block-"+to_string(blocks)+", please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    dofname=ejson.at("dofs").at(i);
                    if(dofname.size()<1){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid or empty in 'dofs' of elmt block-"+to_string(blocks)+", please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    if(!t_dofhandler.isValidDofName(dofname)){
                        MessagePrinter::printErrorTxt("dof-"+to_string(i+1)+"'s name is invalid in 'dofs' of elmt block-"+to_string(blocks)+", it must be one of the names in your 'dofs' block, please check your input file");
                        MessagePrinter::exitAsFem();
                    }
                    elmtBlock.m_DofNames.push_back(dofname);
                    elmtBlock.m_DofIDs.push_back(t_dofhandler.getDofIDViaName(dofname));
                }
                HasDofs=true;
            }
        }
        else{
            MessagePrinter::printErrorTxt("can\'t find dofs name in element block-"+to_string(blocks)+", please check your input file");
            MessagePrinter::exitAsFem();
        }//end-of-'dofs'-reading

        if(ejson.contains("domain")){
            if(ejson.at("domain").size()<1){
                MessagePrinter::printErrorTxt("invalid or empty domain name in '"+elmtBlock.m_ElmtBlockName+"' of your [elmts] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            string domain;
            for(int i=0;i<static_cast<int>(ejson.at("domain").size());i++){
                if(!ejson.at("domain").at(i).is_string()){
                    MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of '"+elmtBlock.m_ElmtBlockName+"' in your 'elmts' subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                domain=ejson.at("domain").at(i);
                if(domain.size()<1){
                    MessagePrinter::printErrorTxt("domain name is invalid or empty in 'domain' of '"+elmtBlock.m_ElmtBlockName+"' in your 'elmts' subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                if(!t_fecell.isFECellBulkElmtPhyNameValid(domain)){
                    MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of '"+elmtBlock.m_ElmtBlockName+"' in your 'elmts' subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                elmtBlock.m_DomainNameList.push_back(domain);
            }
        }
        else{
            elmtBlock.m_DomainNameList.clear();
            elmtBlock.m_DomainNameList.push_back("alldomain");
        } //end-of-'domain'-reading

        if(ejson.contains("material")){
            if(!ejson.at("material").contains("type")){
                MessagePrinter::printErrorTxt("can\'t find material type in element block-"+to_string(blocks)+", please check your input file");
                MessagePrinter::exitAsFem();
            }
            elmtBlock.m_MateTypeName=ejson.at("material").at("type");

            elmtBlock.m_MateType=MateType::NULLMATE;
            if(elmtBlock.m_MateTypeName=="constpoisson"){
                elmtBlock.m_MateType=MateType::CONSTPOISSONMATE;
            }
            else if(elmtBlock.m_MateTypeName=="poisson1d-benchmark"){
                elmtBlock.m_MateType=MateType::POISSON1DBENCHMARKMATE;
            }
            else if(elmtBlock.m_MateTypeName=="poisson2d-benchmark"){
                elmtBlock.m_MateType=MateType::POISSON2DBENCHMARKMATE;
            }
            else if(elmtBlock.m_MateTypeName=="nonlinear-poisson2d"){
                elmtBlock.m_MateType=MateType::NONLINEARPOISSON2DMATE;
            }
            else if(elmtBlock.m_MateTypeName=="nonlinear-poisson3d"){
                elmtBlock.m_MateType=MateType::NONLINEARPOISSON3DMATE;
            }
            else if(elmtBlock.m_MateTypeName=="constdiffusion"){
                elmtBlock.m_MateType=MateType::CONSTDIFFUSIONNMATE;
            }
            else if(elmtBlock.m_MateTypeName=="nonlinear-diffusion2d"){
                elmtBlock.m_MateType=MateType::NONLINEARDIFFUSION2DMATE;
            }
            else if(elmtBlock.m_MateTypeName=="constthermal"){
                elmtBlock.m_MateType=MateType::CONSTTHERMALMATE;
            }
            // for mechanics
            else if(elmtBlock.m_MateTypeName=="linearelastic"){
                elmtBlock.m_MateType=MateType::LINEARELASTICMATE;
            }
            else if(elmtBlock.m_MateTypeName=="saintvenant"){
                elmtBlock.m_MateType=MateType::SAINTVENANTMATE;
            }
            else if(elmtBlock.m_MateTypeName=="neohookean"){
                elmtBlock.m_MateType=MateType::NEOHOOKEANMATE;
            }
            else if(elmtBlock.m_MateTypeName=="smallstrainj2plasticity"){
                elmtBlock.m_MateType=MateType::SMALLSTRAINJ2PLASTICITYMATE;
            }
            else if(elmtBlock.m_MateTypeName=="smallstrainexplawj2plasticity"){
                elmtBlock.m_MateType=MateType::SMALLSTRAINEXPLAWJ2PLASTICITYMATE;
            }
            // for cahn-hilliard material
            else if(elmtBlock.m_MateTypeName=="doublewell"){
                elmtBlock.m_MateType=MateType::DOUBLEWELLMATE;
            }
            else if(elmtBlock.m_MateTypeName=="binarymixture"){
                elmtBlock.m_MateType=MateType::BINARYMIXMATE;
            }
            // for kobayashi dendrite material
            else if(elmtBlock.m_MateTypeName=="kobayashimate"){
                elmtBlock.m_MateType=MateType::KOBAYASHIMATE;
            }
            // for wave material
            else if(elmtBlock.m_MateTypeName=="wave"){
                elmtBlock.m_MateType=MateType::WAVEMATE;
            }
            // for coupled material
            else if(elmtBlock.m_MateTypeName=="smallstraincahnhilliard"){
                elmtBlock.m_MateType=MateType::SMALLSTRAINCAHNHILLIARDMATE;
            }
            else if(elmtBlock.m_MateTypeName=="smallstraindiffusion"){
                elmtBlock.m_MateType=MateType::SMALLSTRAINDIFFUSIONMATE;
            }
            else if(elmtBlock.m_MateTypeName=="smallstraindiffusionj2plasticity"){
                elmtBlock.m_MateType=MateType::SMALLSTRAINDIFFUSIONJ2MATE;
            }
            else if(elmtBlock.m_MateTypeName=="linearelasticfracture"){
                elmtBlock.m_MateType=MateType::LINEARELASTICFRACMATE;
            }
            else if(elmtBlock.m_MateTypeName=="neohookeanpffracture"){
                elmtBlock.m_MateType=MateType::NEOHOOKEANPFFRACTUREMATE;
            }
            else if(elmtBlock.m_MateTypeName=="miehefracture"){
                elmtBlock.m_MateType=MateType::MIEHEFRACTUREMATE;
            }
            else if(elmtBlock.m_MateTypeName=="diffusionacfracture"){
                elmtBlock.m_MateType=MateType::DIFFUSIONACFRACTUREMATE;
            }
            // for user-defined-materials(UMAT)
            else if(elmtBlock.m_MateTypeName=="user1"){
                elmtBlock.m_MateType=MateType::USER1MATE;
            }
            else if(elmtBlock.m_MateTypeName=="user2"){
                elmtBlock.m_MateType=MateType::USER2MATE;
            }
            else if(elmtBlock.m_MateTypeName=="user3"){
                elmtBlock.m_MateType=MateType::USER3MATE;
            }
            else if(elmtBlock.m_MateTypeName=="user4"){
                elmtBlock.m_MateType=MateType::USER4MATE;
            }
            else if(elmtBlock.m_MateTypeName=="user5"){
                elmtBlock.m_MateType=MateType::USER5MATE;
            }
            else if(elmtBlock.m_MateTypeName=="user6"){
                elmtBlock.m_MateType=MateType::USER6MATE;
            }
            else if(elmtBlock.m_MateTypeName=="user7"){
                elmtBlock.m_MateType=MateType::USER7MATE;
            }
            else if(elmtBlock.m_MateTypeName=="user8"){
                elmtBlock.m_MateType=MateType::USER8MATE;
            }
            else if(elmtBlock.m_MateTypeName=="user9"){
                elmtBlock.m_MateType=MateType::USER9MATE;
            }
            else if(elmtBlock.m_MateTypeName=="user10"){
                elmtBlock.m_MateType=MateType::USER10MATE;
            }
            else{
                MessagePrinter::printErrorTxt("unsupported material type name in "+elmtBlock.m_ElmtBlockName," please check your input file");
                MessagePrinter::exitAsFem();
            }
            // we will keep the parameter's json into m_json_params, this will be used later
            // for material property calculation (UMAT)
            if(!ejson.at("material").contains("parameters")){
                MessagePrinter::printWarningTxt("can\'t find parameters in your 'material' of element block-"+to_string(blocks)+", then no parameters will be used");
                elmtBlock.m_JsonParams.clear();
            }
            else{
                elmtBlock.m_JsonParams=ejson.at("material").at("parameters");
            }//for-parameters-in-materials
        }
        else{
            MessagePrinter::printWarningTxt("can\'t find 'material' in element block-"+to_string(blocks)+", then no materials will be used");
            elmtBlock.m_JsonParams.clear();
            elmtBlock.m_MateTypeName="null";
            elmtBlock.m_MateType=MateType::NULLMATE;
        } //end-of-'material'-reading

        if(!HasType || !HasDofs){
            MessagePrinter::printErrorTxt("information in '"+elmtBlock.m_ElmtBlockName
                                         +"' of your 'elmts' subblock is not complete, please check your input file");
            return false;
        }

        t_elmtsystem.addElmtBlock2List(elmtBlock);
    }// end-of-iter-of-'elmts'-block

    return true;
}