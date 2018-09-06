#include "dLaxFriedrichs.h"
//-----------------------------//
dVectorND <double> dLaxFriedrichs::Solve(    const dVectorND <double>& u_m_nP,
                                    const dVectorND <double>& X_mplus1_nP,
                                    const dVectorND <double>& X_mminus1_nP,
                                    const dVectorND <double>& Y_m_nplus1P,
                                    const dVectorND <double>& Y_m_nminus1P,
                                    const dVectorND <double>& F_m_nP,
                                    const dVectorND <double>& R_m_nP,
                                    double DeltaTP,
                                    double DeltaXP,
                                    double DeltaYP) {
    return u_m_nP -
           (X_mplus1_nP - X_mminus1_nP) * (DeltaTP / 2.0 / DeltaXP) -
           (Y_m_nplus1P - Y_m_nminus1P) * (DeltaTP / 2.0 / DeltaYP) -
           F_m_nP -
           R_m_nP;
}