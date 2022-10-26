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
//+++ Date   : 2020.06.28
//+++ Purpose: Implement a general message printer for AsFem
//+++          This class is general to provide almost all the
//+++          message print in AsFem
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#pragma once

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

#include "petsc.h"

#include "Utils/MessageColor.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::to_string;

/**
 * this class offers general message output in the command-line
 */
class MessagePrinter{
public:
    /**
     * constructor
     */
    MessagePrinter();

    /**
     * print out the standard text
     * @param str the string to be printed
     * @param color the color for current string printing
     */
    static void printTxt(string str,MessageColor color=MessageColor::WHITE);

    /**
     * print out a short text string
     * @param str the string to be printed
     * @param color the color for current string printing
     */
    static void printShortTxt(string str,MessageColor color=MessageColor::WHITE);

    /**
     * print out a long text string
     * @param str the string to be printed
     * @param color the color for current string printing
     */
    static void printLongTxt(string str,MessageColor color=MessageColor::WHITE);

    /**
     * print out the text as an error message
     * @param str the string to be printed
     * @param flag true->print out the star in red color, otherwise in white color
     */
    static void printErrorTxt(string str,bool flag=true);
    /**
     * print out the text as a warning message
     * @param str the string to be printed
     * @param flag true->print out the star in yellow color, otherwise in white color
     */
    static void printWarningTxt(string str,bool flag=true);

    /**
     * print out a normal text string
     * @param str the string to be printed
     * @param color the color for current string printing
     */
    static void printNormalTxt(string str,MessageColor color=MessageColor::WHITE);

    /**
     * print out a welcome text
     * @param str the string to be printed
     */
    static void printWelcomeTxt(string str);

    /**
     * print out a star line
     * @param color the color for starts, default is white
     */
    static void printStars(MessageColor color=MessageColor::WHITE);

    /**
     * print out a dashed line
     * @param color the color for the dashed line, default is white
     */
    static void printDashLine(MessageColor color=MessageColor::WHITE);

    //***********************************
    //*** for input file read
    //***********************************
    /**
     * print out an error message indicate the line number, this should be used in your 
     * input file reading
     * @param linenumber integer, the current line of the error line in your input file
     */
    static void printErrorInLineNumber(const int &linenumber);

    /**
     * this will stop all the process and exit the asfem program!
     */
    static void exitAsFem();

    /**
     * set the color which will be used in your terminal output
     * @param color the color for terminal output
     */
    static void setColor(const MessageColor &color);

private:
    static const int _nWords=78;/**< the number of characters in one single line output */

    /**
     * this function split the str into several substring once it is beyond the length limit
     * @param str the input string
     */
    vector<string> splitStr2Vec(string str);
    /**
     * this function split the error str into several substring once it is beyond the length limit
     * @param str the input string
     */
    vector<string> splitErrorStr2Vec(string str);
    /**
     * this function split the warning str into several substring once it is beyond the length limit
     * @param str the input string
     */
    vector<string> splitWarningStr2Vec(string str);
    /**
     * this function split the normal str into several substring once it is beyond the length limit
     * @param str the input string
     */
    vector<string> splitNormalStr2Vec(string str);

};