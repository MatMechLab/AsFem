#ifndef ASFEM_WELCOME_H
#define ASFEM_WELCOME_H

#include <iostream>

void Welcome(const double &version,const int &year,const int &mon,const int &day){
    std::cout<<"**************************************************************"<<std::endl;
    std::cout<<"*** Welcome to use AsFem                                   ***"<<std::endl;
    std::cout<<"*** A Simple Finite Element Method program                 ***"<<std::endl;
    std::printf("*** Version: %6.2f  Release @ %4d.%02d.%02d                  ***\n",version,year,mon,day);
    std::printf("*** Author : walkandthinker @ CopyRight %4d               ***\n",year);
    std::cout<<"*** Contact: walkandthinker@gmail.com                      ***"<<std::endl;
    std::cout<<"*** License: GPL-3.0                                       ***"<<std::endl;
    std::cout<<"*** QQ group: 879908352                                    ***"<<std::endl;
    std::cout<<"*** Feel free to use and discuss .:.                       ***"<<std::endl;
    std::cout<<"**************************************************************"<<std::endl;
}

#endif //ASFEM_WELCOME_H
