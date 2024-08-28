#include <iostream>
#include <fstream>
#include <cmath>

#include "lgd.h"

void LGD::calcEnergyLGD(double *energy, 
                     double *pol,
                     double *alpha,
                     double *beta,
                     double *gamma)
    {
        *energy = 1e0/2e0 * (*alpha) * std::pow(*pol, 2) 
                + 1e0/4e0 * (*beta) * std::pow(*pol, 4)
                + 1e0/6e0 * (*gamma) * std::pow(*pol, 6);
    }

void LGD::calcEnergyLGDDer(double *energyDer,
                        double *pol,
                        double *alpha,
                        double *beta,
                        double *gamma)
    {
        *energyDer = (*alpha) * (*pol)
                 +   (*beta) * std::pow(*pol, 3)
                 +  (*gamma) * std::pow(*pol, 5);
    }

void LGD::calcEnergyLGDDer2(double *energyDer2,
                         double *pol,
                         double *alpha,
                         double *beta,
                         double *gamma)
    {
        *energyDer2 = (*alpha) 
                    + 3e0 * (*beta) * std::pow(*pol, 2)
                    + 5e0 * (*gamma) * std::pow(*pol, 4);
    }

void LGD::calcDFDP(double *DFDP,
               double *pol,
               double *alpha,
               double *beta,
               double *gamma,
               double *E)
    {
        double eps0 = 8.85E-12;
        *DFDP = (*alpha) * (*pol)
              + (*beta) * std::pow(*pol, 3)
              + (*gamma) * std::pow(*pol, 5) - (*E);
        //std::cout << *DFDP << std::endl;
    }

void LGD::LKequation(double *P_new,
                     double *P_old,
                     double *DeltaT,
                     double *gammaVis,
                     double *DFDP)
    {
        *P_new = *P_old - (*DeltaT) * (*gammaVis) * (*DFDP);
        //std::cout << - (*DeltaT) * (*gammaVis) * (*DFDP) << std::endl;
    }
