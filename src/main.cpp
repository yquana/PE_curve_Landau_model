#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ios>

#include "json.hpp"
#include "lgd.h"


int main(int argc, char *argv[])
{
    nlohmann::json data;
    std::fstream infile;
    std::cout << " LINE 15 " << std::endl;
    infile.open(argv[1], std::ios::in);
    data = nlohmann::json::parse(infile);

    double alpha = data["alpha"];
    double beta = data["beta"];
    double gamma = data["gamma"];
    double polInit = data["polInit"];
    double gammaVis = data["gammaVis"];
    double Einit = data["EInit"];
    double EMax  = data["EMax"];
    double EMin  = data["EMin"];
    double EStep = data["EStep"];
    int   nsteps = data["nsteps"];
    double deltaT = data["deltaT"];

    int iStep = 0;
    double pol = polInit;
    double polOld = polInit;
    double polNew = polInit;
    double EField = 0e0;
    double sign = 1E0;
    double energyLGD = 0e0;
    double energyLGDDer = 0e0;
    double energyLGDDer2 = 0e0;
    double energyTotal = 0e0;
    double DFDP = 0e0;

    double eps0 = 8.85E-12;

    if (beta*beta < 4e0 * alpha * gamma)
    {
        std::cout << "The system is not ferroelectric" << std::endl;
    }
    else
    {
        double P2[4];
        P2[0] = (-beta - std::sqrt(beta*beta - 4e0 * alpha * gamma)) / (2e0 * gamma);
        P2[1] = (-beta + std::sqrt(beta*beta - 4e0 * alpha * gamma)) / (2e0 * gamma);
        if (9e0 * beta * beta > 16e0 * alpha * gamma)
        {
            P2[2] = (-3e0 * beta + std::sqrt(9e0 * beta * beta - 16e0 * alpha * gamma)) / (8e0 * gamma);
            P2[3] = (-3e0 * beta - std::sqrt(9e0 * beta * beta - 16e0 * alpha * gamma)) / (8e0 * gamma);
            std::cout << "Pr = " << std::sqrt(P2[0]) << " " << std::sqrt(P2[1]) << std::endl;
            std::cout << "Pc = " << std::sqrt(P2[2]) << " " << std::sqrt(P2[3]) << std::endl;
            std::cout << "E_c = " << alpha * std::sqrt(P2[2]) + beta * std::pow(std::sqrt(P2[2]), 3e0) + gamma * std::pow(std::sqrt(P2[2]), 5e0) << std::endl;
            std::cout << "E_c = " << alpha * std::sqrt(P2[3]) + beta * std::pow(std::sqrt(P2[3]), 3e0) + gamma * std::pow(std::sqrt(P2[3]), 5e0) << std::endl;
        }
        else
        {
            std::cout << "P2 = " << P2[0] << " " << P2[1] << std::endl;
        }
    }

    std::fstream outfile;
    std::string ofname;
    ofname = std::string(argv[1]) + "output.dat";
    outfile.open(ofname.c_str(), std::ios::out);
    outfile << "nstep EField pol DFDP" << std::endl;
    while(iStep < nsteps)
    {
        if(EField > EMax || EField < EMin)
        {
            sign *= -1E0;
        }
        for(int l = 0; l < 10; l++)
        {
            LGD::calcDFDP(&DFDP, &polOld, &alpha, &beta, &gamma, &EField);
            LGD::LKequation(&polNew, &polOld, &deltaT, &gammaVis, &DFDP);
            polOld = polNew;
        }
        EField += sign * EStep;
        outfile << iStep << " " << EField << " " << polNew << " " << DFDP <<  std::endl;
        iStep++;
    }
    outfile.close();
}
