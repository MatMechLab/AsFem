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
//+++ Date   : 2022.06.05
//+++ Purpose: implement the 1d lobatto gauss point generator
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "FE/QPoint1DLobattoGenerator.h"

void QPoint1DLobattoGenerator::generateQPoints(const int &t_order,const MeshType &t_meshtype,int &t_ngp,vector<double> &t_qpoints){
    if(t_meshtype==MeshType::EDGE2){}// for 1d case, we don't need to consider the mesh type
    // the numbers are generated from:
    //    https://keisan.casio.com/exec/system/1280801905
    // here order=2*n-3, which means it needs at least two points
    switch (t_order)
    {
    case 0:
    case 1:{
        t_ngp=2;
        t_qpoints.resize(t_ngp*2,0.0);
        t_qpoints[(1-1)*2+0]= 1.0;
        t_qpoints[(1-1)*2+1]=-1.0;

        t_qpoints[(2-1)*2+0]= 1.0;
        t_qpoints[(2-1)*2+1]= 1.0;
        break;
        }
    case 2:
    case 3:{
        t_ngp=3;
        t_qpoints.resize(t_ngp*2,0.0);
        t_qpoints[(1-1)*2+0]= 1.0/3.0;
        t_qpoints[(1-1)*2+1]=-1.0;

        t_qpoints[(2-1)*2+0]= 4.0/3.0;
        t_qpoints[(2-1)*2+1]= 0.0;

        t_qpoints[(3-1)*2+0]= 1.0/3.0;
        t_qpoints[(3-1)*2+1]= 1.0;
        break;
        }
    case 4:
    case 5:{
        t_ngp=4;
        t_qpoints.resize(t_ngp*2,0.0);
        t_qpoints[(1-1)*2+0]= 1.0/6.0;
        t_qpoints[(1-1)*2+1]=-1.0;

        t_qpoints[(2-1)*2+0]= 5.0/6.0;
        t_qpoints[(2-1)*2+1]=-sqrt(1.0/5.0);

        t_qpoints[(3-1)*2+0]= 5.0/6.0;
        t_qpoints[(3-1)*2+1]= sqrt(1.0/5.0);

        t_qpoints[(4-1)*2+0]= 1.0/6.0;
        t_qpoints[(4-1)*2+1]= 1.0;
        break;
        }
    case 6:
    case 7:{
        t_ngp=5;
        t_qpoints.resize(t_ngp*2,0.0);
        t_qpoints[(1-1)*2+0]= 0.1;
        t_qpoints[(1-1)*2+1]=-1.0;

        t_qpoints[(2-1)*2+0]= 49.0/90.0;
        t_qpoints[(2-1)*2+1]=-sqrt(3.0/7.0);

        t_qpoints[(3-1)*2+0]= 32.0/45.0;
        t_qpoints[(3-1)*2+1]= 0.0;

        t_qpoints[(4-1)*2+0]= 49.0/90.0;
        t_qpoints[(4-1)*2+1]= sqrt(3.0/7.0);

        t_qpoints[(5-1)*2+0]= 0.1;
        t_qpoints[(5-1)*2+1]= 1.0;

        break;
        }
    case 8:
    case 9:{
        t_ngp=6;
        t_qpoints.resize(t_ngp*2,0.0);

        t_qpoints[(1-1)*2+0]= 0.06666666666666666666667;
        t_qpoints[(1-1)*2+1]=-1.0;

        t_qpoints[(2-1)*2+0]= 0.3784749562978469803166;
        t_qpoints[(2-1)*2+1]=-0.765055323929464692851;

        t_qpoints[(3-1)*2+0]= 0.5548583770354863530167;
        t_qpoints[(3-1)*2+1]=-0.2852315164806450963142;

        t_qpoints[(4-1)*2+0]= 0.5548583770354863530167;
        t_qpoints[(4-1)*2+1]= 0.2852315164806450963142;

        t_qpoints[(5-1)*2+0]= 0.3784749562978469803166;
        t_qpoints[(5-1)*2+1]= 0.765055323929464692851;

        t_qpoints[(6-1)*2+0]= 0.06666666666666666666667;
        t_qpoints[(6-1)*2+1]= 1.0;
        
        break;
        }
    case 10:
    case 11:{
        t_ngp=7;
        t_qpoints.resize(t_ngp*2,0.0);

        t_qpoints[(1-1)*2+0]= 0.04761904761904761904762;
        t_qpoints[(1-1)*2+1]=-1.0;

        t_qpoints[(2-1)*2+0]= 0.276826047361565948011;
        t_qpoints[(2-1)*2+1]=-0.830223896278566929872;

        t_qpoints[(3-1)*2+0]= 0.4317453812098626234179;
        t_qpoints[(3-1)*2+1]=-0.4688487934707142138038;

        t_qpoints[(4-1)*2+0]= 0.487619047619047619048;
        t_qpoints[(4-1)*2+1]= 0.0;

        t_qpoints[(5-1)*2+0]= 0.431745381209862623418;
        t_qpoints[(5-1)*2+1]= 0.468848793470714213804;

        t_qpoints[(6-1)*2+0]= 0.2768260473615659480107;
        t_qpoints[(6-1)*2+1]= 0.830223896278566929872;

        t_qpoints[(7-1)*2+0]= 0.04761904761904761904762;
        t_qpoints[(7-1)*2+1]= 1.0;
        
        break;
        }
    case 12:
    case 13:{
        t_ngp=8;
        t_qpoints.resize(t_ngp*2,0.0);

        t_qpoints[(1-1)*2+0]= 0.03571428571428571428571;
        t_qpoints[(1-1)*2+1]=-1.0;

        t_qpoints[(2-1)*2+0]= 0.210704227143506039383;
        t_qpoints[(2-1)*2+1]=-0.8717401485096066153375;

        t_qpoints[(3-1)*2+0]= 0.3411226924835043647642;
        t_qpoints[(3-1)*2+1]=-0.5917001814331423021445;

        t_qpoints[(4-1)*2+0]= 0.4124587946587038815671;
        t_qpoints[(4-1)*2+1]=-0.2092992179024788687687;

        t_qpoints[(5-1)*2+0]= 0.412458794658703881567;
        t_qpoints[(5-1)*2+1]= 0.2092992179024788687687;

        t_qpoints[(6-1)*2+0]= 0.341122692483504364764;
        t_qpoints[(6-1)*2+1]= 0.5917001814331423021445;

        t_qpoints[(7-1)*2+0]= 0.210704227143506039383;
        t_qpoints[(7-1)*2+1]= 0.8717401485096066153375;

        t_qpoints[(8-1)*2+0]= 0.03571428571428571428571;
        t_qpoints[(8-1)*2+1]= 1.0;
        
        break;
        }
    case 14:
    case 15:{
        t_ngp=9;
        t_qpoints.resize(t_ngp*2,0.0);

        t_qpoints[(1-1)*2+0]= 0.02777777777777777777778;
        t_qpoints[(1-1)*2+1]=-1.0;

        t_qpoints[(2-1)*2+0]= 0.1654953615608055250463;
        t_qpoints[(2-1)*2+1]=-0.8997579954114601573124;

        t_qpoints[(3-1)*2+0]= 0.274538712500161735281;
        t_qpoints[(3-1)*2+1]=-0.6771862795107377534459;

        t_qpoints[(4-1)*2+0]= 0.3464285109730463451151;
        t_qpoints[(4-1)*2+1]=-0.3631174638261781587108;

        t_qpoints[(5-1)*2+0]= 0.3715192743764172335601;
        t_qpoints[(5-1)*2+1]= 0.0;

        t_qpoints[(6-1)*2+0]= 0.3464285109730463451151;
        t_qpoints[(6-1)*2+1]= 0.3631174638261781587108;

        t_qpoints[(7-1)*2+0]= 0.2745387125001617352807;
        t_qpoints[(7-1)*2+1]= 0.6771862795107377534459;

        t_qpoints[(8-1)*2+0]= 0.165495361560805525046;
        t_qpoints[(8-1)*2+1]= 0.8997579954114601573124;

        t_qpoints[(9-1)*2+0]= 0.02777777777777777777778;
        t_qpoints[(9-1)*2+1]= 1.0;
        
        break;
        }
    case 16:
    case 17:{
        t_ngp=10;
        t_qpoints.resize(t_ngp*2,0.0);

        t_qpoints[(1-1)*2+0]= 0.02222222222222222222222;
        t_qpoints[(1-1)*2+1]=-1.0;

        t_qpoints[(2-1)*2+0]= 0.1333059908510701111262;
        t_qpoints[(2-1)*2+1]=-0.9195339081664588138289;

        t_qpoints[(3-1)*2+0]= 0.2248893420631264521195;
        t_qpoints[(3-1)*2+1]=-0.7387738651055050750031;

        t_qpoints[(4-1)*2+0]= 0.2920426836796837578756;
        t_qpoints[(4-1)*2+1]=-0.4779249498104444956612;

        t_qpoints[(5-1)*2+0]= 0.3275397611838974566565;
        t_qpoints[(5-1)*2+1]=-0.1652789576663870246262;

        t_qpoints[(6-1)*2+0]= 0.3275397611838974566565;
        t_qpoints[(6-1)*2+1]= 0.1652789576663870246262;

        t_qpoints[(7-1)*2+0]= 0.292042683679683757876;
        t_qpoints[(7-1)*2+1]= 0.4779249498104444956612;

        t_qpoints[(8-1)*2+0]= 0.224889342063126452119;
        t_qpoints[(8-1)*2+1]= 0.7387738651055050750031;

        t_qpoints[(9-1)*2+0]= 0.133305990851070111126;
        t_qpoints[(9-1)*2+1]= 0.9195339081664588138289;

        t_qpoints[(10-1)*2+0]= 0.02222222222222222222222;
        t_qpoints[(10-1)*2+1]= 1.0;
        
        break;
        }
    default:{
        MessagePrinter::printErrorTxt("order="+to_string(t_order)+" is not supported in QPoint1DLobattoGenerator");
        MessagePrinter::exitAsFem();
        break;
        }
    }
}