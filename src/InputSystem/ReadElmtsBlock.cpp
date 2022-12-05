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
//+++ Date   : 2022.05.12
//+++ Purpose: read the elements block from json
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "InputSystem/InputSystem.h"
#include "ElmtSystem/ElmtBlock.h"

bool InputSystem::readElmtsBlock(nlohmann::json &t_json,const Mesh &t_mesh,const DofHandler &t_dofhandler,ElmtSystem &t_elmtsystem){
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
        
        elmtBlock.m_elmtblock_index=blocks;
        elmtBlock.m_elmt_blockname=it.key();
        nlohmann::json ejson=t_json.at(elmtBlock.m_elmt_blockname);//it.value();// now ejson gose into 'elmt-x'
        if(ejson.contains("type")){
            if(!ejson.at("type").is_string()){
                MessagePrinter::printErrorTxt("invalid type name in element block-"+to_string(blocks)+", please check your input file");
                MessagePrinter::exitAsFem();
            }
            elmtBlock.m_elmt_typename=ejson.at("type");
            elmtBlock.m_elmttype=ElmtType::NULLELMT;
            if(elmtBlock.m_elmt_typename=="poisson"){
                elmtBlock.m_elmttype=ElmtType::POISSONELMT;
            }
            else if(elmtBlock.m_elmt_typename=="diffusion"){
                elmtBlock.m_elmttype=ElmtType::DIFFUSIONELMT;
            }
            else if(elmtBlock.m_elmt_typename=="thermal"){
                elmtBlock.m_elmttype=ElmtType::THERMALELMT;
            }
            else if(elmtBlock.m_elmt_typename=="mechanics"){
                elmtBlock.m_elmttype=ElmtType::MECHANICSELMT;
            }
            else if(elmtBlock.m_elmt_typename=="dynamics"){
                elmtBlock.m_elmttype=ElmtType::DYNAMICMECHANICSELMT;
            }
            else if(elmtBlock.m_elmt_typename=="allencahn"){
                elmtBlock.m_elmttype=ElmtType::ALLENCAHNELMT;
            }
            else if(elmtBlock.m_elmt_typename=="cahnhilliard"){
                elmtBlock.m_elmttype=ElmtType::CAHNHILLIARDELMT;
            }
            else if(elmtBlock.m_elmt_typename=="wave"){
                elmtBlock.m_elmttype=ElmtType::WAVEELMT;
            }
            else if(elmtBlock.m_elmt_typename=="kobayashi"){
                elmtBlock.m_elmttype=ElmtType::KOBAYASHIELMT;
            }
            else if(elmtBlock.m_elmt_typename=="laplace"){
                elmtBlock.m_elmttype=ElmtType::LAPLACEELMT;
            }
            else if(elmtBlock.m_elmt_typename=="scalarbodysource"){
                elmtBlock.m_elmttype=ElmtType::SCALARBODYSOURCEELMT;
            }
            // for coupled cases
            else if(elmtBlock.m_elmt_typename=="stresscahnhilliard"){
                elmtBlock.m_elmttype=ElmtType::STRESSCAHNHILLIARDELMT;
            }
            else if(elmtBlock.m_elmt_typename=="stressdiffusion"){
                elmtBlock.m_elmttype=ElmtType::STRESSDIFFUSIONELMT;
            }
            else if(elmtBlock.m_elmt_typename=="allencahnfracture"){
                elmtBlock.m_elmttype=ElmtType::ALLENCAHNFRACTUREELMT;
            }
            else if(elmtBlock.m_elmt_typename=="miehefracture"){
                elmtBlock.m_elmttype=ElmtType::MIEHEFRACTUREELMT;
            }
            else if(elmtBlock.m_elmt_typename=="diffusionacfracture"){
                elmtBlock.m_elmttype=ElmtType::DIFFUSIONACFRACTUREELMT;
            }
            // for user-defined-element
            else if(elmtBlock.m_elmt_typename=="user1"){
                elmtBlock.m_elmttype=ElmtType::USER1ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user2"){
                elmtBlock.m_elmttype=ElmtType::USER2ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user3"){
                elmtBlock.m_elmttype=ElmtType::USER3ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user4"){
                elmtBlock.m_elmttype=ElmtType::USER4ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user5"){
                elmtBlock.m_elmttype=ElmtType::USER5ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user6"){
                elmtBlock.m_elmttype=ElmtType::USER6ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user7"){
                elmtBlock.m_elmttype=ElmtType::USER7ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user8"){
                elmtBlock.m_elmttype=ElmtType::USER8ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user9"){
                elmtBlock.m_elmttype=ElmtType::USER9ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user10"){
                elmtBlock.m_elmttype=ElmtType::USER10ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user11"){
                elmtBlock.m_elmttype=ElmtType::USER11ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user12"){
                elmtBlock.m_elmttype=ElmtType::USER12ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user13"){
                elmtBlock.m_elmttype=ElmtType::USER13ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user14"){
                elmtBlock.m_elmttype=ElmtType::USER14ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user15"){
                elmtBlock.m_elmttype=ElmtType::USER15ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user16"){
                elmtBlock.m_elmttype=ElmtType::USER16ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user17"){
                elmtBlock.m_elmttype=ElmtType::USER17ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user18"){
                elmtBlock.m_elmttype=ElmtType::USER18ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user19"){
                elmtBlock.m_elmttype=ElmtType::USER19ELMT;
            }
            else if(elmtBlock.m_elmt_typename=="user20"){
                elmtBlock.m_elmttype=ElmtType::USER20ELMT;
            }
            else{
                MessagePrinter::printErrorTxt("unsupported element type in "+elmtBlock.m_elmt_blockname+", please check your input file");
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
                    elmtBlock.m_dof_names.push_back(dofname);
                    elmtBlock.m_dof_ids.push_back(t_dofhandler.getDofIDViaName(dofname));
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
                MessagePrinter::printErrorTxt("invalid or empty domain name in ["+elmtBlock.m_elmt_blockname+"] of your [elmts] subblock, please check your input file");
                MessagePrinter::exitAsFem();
            }
            string domain;
            for(int i=0;i<static_cast<int>(ejson.at("domain").size());i++){
                if(!ejson.at("domain").at(i).is_string()){
                    MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of ["+elmtBlock.m_elmt_blockname+"] in your [elmts] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                domain=ejson.at("domain").at(i);
                if(domain.size()<1){
                    MessagePrinter::printErrorTxt("domain name is invalid or empty in 'domain' of ["+elmtBlock.m_elmt_blockname+"] in your [elmts] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                if(!t_mesh.isBulkElmtPhyNameValid(domain)){
                   MessagePrinter::printErrorTxt("domain name is invalid in 'domain' of ["+elmtBlock.m_elmt_blockname+"] in your [elmts] subblock, please check your input file");
                    MessagePrinter::exitAsFem();
                }
                elmtBlock.m_domain_namelist.push_back(domain);
            }
        }
        else{
            elmtBlock.m_domain_namelist.clear();
            elmtBlock.m_domain_namelist.push_back("alldomain");
        } //end-of-'domain'-reading

        if(ejson.contains("material")){
            if(!ejson.at("material").contains("type")){
                MessagePrinter::printErrorTxt("can\'t find material type in element block-"+to_string(blocks)+", please check your input file");
                MessagePrinter::exitAsFem();
            }
            elmtBlock.m_mate_typename=ejson.at("material").at("type");

            elmtBlock.m_matetype=MateType::NULLMATE;
            if(elmtBlock.m_mate_typename=="constpoisson"){
                elmtBlock.m_matetype=MateType::CONSTPOISSONMATE;
            }
            else if(elmtBlock.m_mate_typename=="poisson1d-benchmark"){
                elmtBlock.m_matetype=MateType::POISSON1DBENCHMARKMATE;
            }
            else if(elmtBlock.m_mate_typename=="poisson2d-benchmark"){
                elmtBlock.m_matetype=MateType::POISSON2DBENCHMARKMATE;
            }
            else if(elmtBlock.m_mate_typename=="nonlinear-poisson2d"){
                elmtBlock.m_matetype=MateType::NONLINEARPOISSON2DMATE;
            }
            else if(elmtBlock.m_mate_typename=="nonlinear-poisson3d"){
                elmtBlock.m_matetype=MateType::NONLINEARPOISSON3DMATE;
            }
            else if(elmtBlock.m_mate_typename=="constdiffusion"){
                elmtBlock.m_matetype=MateType::CONSTDIFFUSIONNMATE;
            }
            else if(elmtBlock.m_mate_typename=="nonlinear-diffusion2d"){
                elmtBlock.m_matetype=MateType::NONLINEARDIFFUSION2DMATE;
            }
            else if(elmtBlock.m_mate_typename=="constthermal"){
                elmtBlock.m_matetype=MateType::CONSTTHERMALMATE;
            }
            // for mechanics
            else if(elmtBlock.m_mate_typename=="linearelastic"){
                elmtBlock.m_matetype=MateType::LINEARELASTICMATE;
            }
            else if(elmtBlock.m_mate_typename=="saintvenant"){
                elmtBlock.m_matetype=MateType::SAINTVENANTMATE;
            }
            else if(elmtBlock.m_mate_typename=="neohookean"){
                elmtBlock.m_matetype=MateType::NEOHOOKEANMATE;
            }
            else if(elmtBlock.m_mate_typename=="smallstrainj2plasticity"){
                elmtBlock.m_matetype=MateType::SMALLSTRAINJ2PLASTICITYMATE;
            }
            else if(elmtBlock.m_mate_typename=="smallstrainexplawj2plasticity"){
                elmtBlock.m_matetype=MateType::SMALLSTRAINEXPLAWJ2PLASTICITYMATE;
            }
            // for cahn-hilliard material
            else if(elmtBlock.m_mate_typename=="doublewell"){
                elmtBlock.m_matetype=MateType::DOUBLEWELLMATE;
            }
            else if(elmtBlock.m_mate_typename=="binarymixture"){
                elmtBlock.m_matetype=MateType::BINARYMIXMATE;
            }
            // for kobayashi dendrite material
            else if(elmtBlock.m_mate_typename=="kobayashimate"){
                elmtBlock.m_matetype=MateType::KOBAYASHIMATE;
            }
            // for wave material
            else if(elmtBlock.m_mate_typename=="wave"){
                elmtBlock.m_matetype=MateType::WAVEMATE;
            }
            // for coupled material
            else if(elmtBlock.m_mate_typename=="smallstraincahnhilliard"){
                elmtBlock.m_matetype=MateType::SMALLSTRAINCAHNHILLIARDMATE;
            }
            else if(elmtBlock.m_mate_typename=="smallstraindiffusion"){
                elmtBlock.m_matetype=MateType::SMALLSTRAINDIFFUSIONMATE;
            }
            else if(elmtBlock.m_mate_typename=="smallstraindiffusionj2plasticity"){
                elmtBlock.m_matetype=MateType::SMALLSTRAINDIFFUSIONJ2MATE;
            }
            else if(elmtBlock.m_mate_typename=="linearelasticfracture"){
                elmtBlock.m_matetype=MateType::LINEARELASTICFRACMATE;
            }
            else if(elmtBlock.m_mate_typename=="neohookeanpffracture"){
                elmtBlock.m_matetype=MateType::NEOHOOKEANPFFRACTUREMATE;
            }
            else if(elmtBlock.m_mate_typename=="miehefracture"){
                elmtBlock.m_matetype=MateType::MIEHEFRACTUREMATE;
            }
            else if(elmtBlock.m_mate_typename=="diffusionacfracture"){
                elmtBlock.m_matetype=MateType::DIFFUSIONACFRACTUREMATE;
            }
            // for user-defined-materials(UMAT)
            else if(elmtBlock.m_mate_typename=="user1"){
                elmtBlock.m_matetype=MateType::USER1MATE;
            }
            else if(elmtBlock.m_mate_typename=="user2"){
                elmtBlock.m_matetype=MateType::USER2MATE;
            }
            else if(elmtBlock.m_mate_typename=="user3"){
                elmtBlock.m_matetype=MateType::USER3MATE;
            }
            else if(elmtBlock.m_mate_typename=="user4"){
                elmtBlock.m_matetype=MateType::USER4MATE;
            }
            else if(elmtBlock.m_mate_typename=="user5"){
                elmtBlock.m_matetype=MateType::USER5MATE;
            }
            else if(elmtBlock.m_mate_typename=="user6"){
                elmtBlock.m_matetype=MateType::USER6MATE;
            }
            else if(elmtBlock.m_mate_typename=="user7"){
                elmtBlock.m_matetype=MateType::USER7MATE;
            }
            else if(elmtBlock.m_mate_typename=="user8"){
                elmtBlock.m_matetype=MateType::USER8MATE;
            }
            else if(elmtBlock.m_mate_typename=="user9"){
                elmtBlock.m_matetype=MateType::USER9MATE;
            }
            else if(elmtBlock.m_mate_typename=="user10"){
                elmtBlock.m_matetype=MateType::USER10MATE;
            }
            else if(elmtBlock.m_mate_typename=="user11"){
                elmtBlock.m_matetype=MateType::USER11MATE;
            }
            else if(elmtBlock.m_mate_typename=="user12"){
                elmtBlock.m_matetype=MateType::USER12MATE;
            }
            else if(elmtBlock.m_mate_typename=="user13"){
                elmtBlock.m_matetype=MateType::USER13MATE;
            }
            else if(elmtBlock.m_mate_typename=="user14"){
                elmtBlock.m_matetype=MateType::USER14MATE;
            }
            else if(elmtBlock.m_mate_typename=="user15"){
                elmtBlock.m_matetype=MateType::USER15MATE;
            }
            else if(elmtBlock.m_mate_typename=="user16"){
                elmtBlock.m_matetype=MateType::USER16MATE;
            }
            else if(elmtBlock.m_mate_typename=="user17"){
                elmtBlock.m_matetype=MateType::USER17MATE;
            }
            else if(elmtBlock.m_mate_typename=="user18"){
                elmtBlock.m_matetype=MateType::USER18MATE;
            }
            else if(elmtBlock.m_mate_typename=="user19"){
                elmtBlock.m_matetype=MateType::USER19MATE;
            }
            else if(elmtBlock.m_mate_typename=="user20"){
                elmtBlock.m_matetype=MateType::USER20MATE;
            }
            else{
                MessagePrinter::printErrorTxt("unsupported material type name in "+elmtBlock.m_elmt_blockname," please check your input file");
                MessagePrinter::exitAsFem();
            }
            // we will keep the parameter's json into m_json_params, this will be used later
            // for material property calculation (UMAT)
            if(!ejson.at("material").contains("parameters")){
                MessagePrinter::printWarningTxt("can\'t find parameters in your 'material' of element block-"+to_string(blocks)+", then no parameters will be used");
                elmtBlock.m_json_params.clear();
            }
            else{
                elmtBlock.m_json_params=ejson.at("material").at("parameters");
            }//for-parameters-in-materials
        }
        else{
            MessagePrinter::printWarningTxt("can\'t find 'material' in element block-"+to_string(blocks)+", then no materials will be used");
            elmtBlock.m_json_params.clear();
            elmtBlock.m_mate_typename="null";
            elmtBlock.m_matetype=MateType::NULLMATE;
        } //end-of-'material'-reading

        if(!HasType || !HasDofs){
            MessagePrinter::printErrorTxt("information in ["+elmtBlock.m_elmt_blockname
                                         +"] of your [elmts] subblock is not complete, please check your input file");
            return false;
        }

        t_elmtsystem.addElmtBlock2List(elmtBlock);
    }// end-of-iter-of-'elmts'-block

    return true;
}