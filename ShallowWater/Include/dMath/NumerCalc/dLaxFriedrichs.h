#ifndef DLAXFRIEDRICHS_H
#define DLAXFRIEDRICHS_H
//-----------------------------//
#include "dVectors.h"
//-----------------------------//
class dLaxFriedrichs {
public:
    dLaxFriedrichs() = default;
    ~dLaxFriedrichs() = default;

    dVectorND <double> Solve(   const dVectorND <double>& u_m_nP,
                                const dVectorND <double>& X_mplus1_nP,
                                const dVectorND <double>& X_mminus1_nP,
                                const dVectorND <double>& Y_m_nplus1P,
                                const dVectorND <double>& Y_m_nminus1P,
                                const dVectorND <double>& F_m_nP,
                                const dVectorND <double>& R_m_nP,
                                double DeltaTP,
                                double DeltaXP,
                                double DeltaYP);
private:
};
//-----------------------------//
#endif
