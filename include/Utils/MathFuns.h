//****************************************************************
//* This file is part of the AsFem framework
//* A Simple Finite Element Method program (AsFem)
//* All rights reserved, Yang Bai/M3 Group @ CopyRight 2022
//* https://github.com/M3Group/AsFem
//* Licensed under GNU GPLv3, please see LICENSE for details
//* https://www.gnu.org/licenses/gpl-3.0.en.html
//****************************************************************
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++ Author : Yang Bai
//+++ Date   : 2021.01.17
//+++ Purpose: Implement some commonly used mathematic functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>

#include "Utils/MessagePrinter.h"

inline double BracketPos(const double &x){
    return 0.5*(x+abs(x));
}
//******************************************
inline double BracketNeg(const double &x){
    return 0.5*(x-abs(x));
}

inline double PicewiseLinearInterpolation(const vector<double> &xs,const vector<double> &ys,const double &xi){
    // this function responsible for the picewise linear interpolation function
    // here the x must be a monotonically non-decreasing vector, its starting point is always 0.0!!!
    if(static_cast<int>(xs.size())<1){
        MessagePrinter::PrintErrorTxt("invalid usage of picewise linear interpolation function, your x vector size is <1!!!");
        MessagePrinter::AsFem_Exit();
    }
    double x1,x2,y1,y2,y,a,b;
    bool FoundPoint=false;
    for(int i=0;i<static_cast<int>(xs.size())-1;i++){
        x1=xs[i];x2=xs[i+1];
        y1=ys[i];y2=ys[i+1];
        if(xi>=x1 && xi<x2){
            FoundPoint=true;
            // the line is y=a*x+b
            a=(y2-y1)/(x2-x1);
            b=y1-a*x1;
            y=a*xi+b;
        }
    }
    if(FoundPoint){
        return y;
    }
    else{
        int i=static_cast<int>(xs.size());
        if(xi>xs[i-1]){
            x1=xs[i-2];x2=xs[i-1];
            y1=ys[i-2];y2=ys[i-1];
            a=(y2-y1)/(x2-x1);
            b=y1-a*x1;
            y=a*xi+b;
            return y;
        }
        else{
            MessagePrinter::PrintErrorTxt("can not found the [x1,x2] period for your xi in the picewise linear interpolation");
            MessagePrinter::AsFem_Exit();
        }
    }
    return 0.0;
}
