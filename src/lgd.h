
namespace LGD
{
    // F = 1/2 * alpha * pol^2 + 1/4 * beta * pol^4 + 1/6 * gamma * pol^6 - E * P
    void calcEnergyLGD(double *energyLGD, 
                    double *pol,
                    double *alpha,
                    double *beta,
                    double *gamma);

    void calcEnergyLGDDer(double *energyLGDDer,
                       double *pol,
                       double *alpha,
                       double *beta,
                       double *gamma);
    
    void calcEnergyLGDDer2(double *energyLGDDer2,
                        double *pol,
                        double *alpha,
                        double *beta,
                        double *gamma);

    void calcDFDP(double *DFDP,
                  double *pol,
                  double *alpha,
                  double *beta,
                  double *gamma,
                  double *E);    

    void LKequation(double *P_new, 
                    double *P_old,
                    double *DeltaT,
                    double *gammaVis,
                    double *DFDP);
}