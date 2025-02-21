#include "power_law.h"

using namespace std;
using namespace H5;


void power_law(int &timestep, int &ynpl, double &ddd){
        if (timestep > tyield) {
            for (int i = 0; i < Ny; i++) {
                // if (y(i)  >= upper_block && y(i)  <= lower_block) {
                if (i == line_fault) {
                    for (int j = 0; j < Nx; j++) {
                        // reducing matrix calls
                        const double arsf_temp = ARSF(i, j);
                        const double brsf_temp = BRSF(i, j);
                        const double lrsf_temp = LRSF(i, j);
                        const double fric_temp = FRIC(i, j);
                        const double eta0_temp = ETA0(i, j);

                        const double sxx_temp = SXX.block(i, j, 2, 2).sum();
                        const double syy_temp = SYY.block(i, j, 2, 2).sum();

                        // SXX, pt are averaged from four surrounding pressure nodes
                        if (antiplane) {
                            SIIB(i, j) = sqrt(.5 * (pow(SXX(i, j), 2) + 
                                                    pow(SYY(i, j), 2)) + 
                                                    pow(SXY(i, j), 2) + 
                                                    pow(SZX(i, j), 2) + 
                                                    pow(SZY(i, j), 2)
                                                );                        
                        } else {
                            SIIB(i, j)    = sqrt(pow(SXY(i, j), 2) + 
                                                (pow(sxx_temp, 2) + 
                                                pow(syy_temp, 2) + 
                                                pow(sxx_temp + syy_temp, 2))/ 32.
                                            );
                        }
                        // Computing "elastic" stress invariant
                        const double siiel = SIIB(i, j) / (ETA(i, j) / (GGG(i, j) * dt + ETA(i, j)));
                        
                        // Compute old viscoplastic slip rate
                        // Compute PEFF
                        const double prB_temp = (pt.block(i, j, 2, 2).sum() - pf.block(i, j, 2, 2).sum()) / 4.;
                        double prB = prB_temp;
                        
                        if (prB < 1e3) {
                            prB = 1e3;
                        }
                        // Compute old power law strain rate
                        double SIIB1 = SIIB(i, j);
                        
                        // Compute slip velocity for current stress invariant and state
                        double EIISLIP = V0 * sinh(max(SIIB1, 0.) / arsf_temp / prB) * exp(-(brsf_temp * OM(i, j) + fric_temp) / arsf_temp) / dx;
                        
                        // Compute new ETAVP
                        double ETAPL        = SIIB1 / 2. / EIISLIP;
                        double ETAVP        = 1. / (1. / eta0_temp + 1. / ETAPL);
                        // Compute new stress invariant
                        double SIIB2        = siiel * ETAVP / (GGG(i, j) * dt + ETAVP);
                        const double DSIIB1 = SIIB2 - SIIB1;
                        
                        // Compute slip velocity for current stress invariant and state
                        double V            = 2 * V0 * sinh(max(SIIB2, 0.) / arsf_temp / prB) * 
                                            exp(-(brsf_temp * OM(i, j) + fric_temp) / arsf_temp);
                        
                        EIISLIP = V / dx / 2.;
                        
                        // Compute new ETAVP
                        ETAPL               = SIIB2 / 2. / EIISLIP;
                        ETAVP               = 1. / (1. / eta0_temp + 1. / ETAPL);
                        // Compute new stress invariant
                        const double DSIIB2 = siiel * ETAVP / (GGG(i, j) * dt + ETAVP) - SIIB2;
                        double SIIB4        = 0.;
                        
                        if ((DSIIB1 >= 0 && DSIIB2 <= 0) || (DSIIB1 <= 0 && DSIIB2 >= 0)) {
                            double DSIIB = 1e9;
                            
                            while(abs(DSIIB) > 1e-3) {
                                SIIB4   = (SIIB1 + SIIB2) / 2.;
                                
                                //Compute slip velocity for current stress invariant and state
                                V       = 2 * V0 * sinh(max((SIIB4), 0.) / arsf_temp / prB) * 
                                        exp( - (brsf_temp * OM(i, j) + fric_temp) / arsf_temp);
                                
                                EIISLIP = V / dx / 2.;
                                
                                // Compute new ETAVP
                                ETAPL   = SIIB4 / 2. / EIISLIP;
                                ETAVP   = 1. / (1. / eta0_temp + 1. / ETAPL);
                                // Compute new stress invariant
                                DSIIB   = siiel * ETAVP / (GGG(i, j) * dt + ETAVP) - SIIB4;
                                if ((DSIIB >= 0 && DSIIB1 >= 0) || (DSIIB <= 0 && DSIIB1 <= 0)) {
                                    SIIB1 = SIIB4;
                                } else {
                                    SIIB2 = SIIB4;
                                }
                            }
                        }
                        
                        if (V * dt / lrsf_temp > 1e-6) {
                            OM5(i, j) = log(V0 / V + (exp(OM0(i, j)) - V0 / V) * exp( - V * dt / lrsf_temp));
                        } else {
                            OM5(i, j) = log(exp(OM0(i, j)) * (1 - V * dt / lrsf_temp) + V0 * dt / lrsf_temp);
                        }
                        
                        // Compute yielding stress
                        const double syield = max(syieldmin, prB_temp * arsf_temp * asinh(V / 2. / V0 * exp((brsf_temp * OM5(i, j) + fric_temp) / arsf_temp)));
                        
                        // Compute visco - plastic viscosity
                        double etapl = eta0_temp * syield / (eta0_temp * V + syield);
                        
                        // Save syield
                        SigmaY(i, j)        = syield;
                        VSLIPB(i, j)        = V;
                        SII_fault(i, j)     = SIIB4;

                        // reduces calls on matrix
                        const double g_temp = 2 * GGG(i, j); // reduces calls on matrix
                        
                        // "/ BETASOLID" -> Timestep criterion, Lapusta et al., 2000; Lapusta and Liu, 2009
                        const double k = g_temp / (pi * (1 - ((3. / BETASOLID - g_temp) / (6. / BETASOLID + g_temp))) * dx);
                        double dTETAmax;
                        if ((pow((k * lrsf_temp / prB - brsf_temp) / arsf_temp - 1, 2) / 4. - k * lrsf_temp / arsf_temp / prB) < 0) {
                            dTETAmax = min(1. - (brsf_temp - arsf_temp) * prB / (k * lrsf_temp), .2);
                        } else {
                            dTETAmax = min(arsf_temp * prB / (k * lrsf_temp - (brsf_temp - arsf_temp) * prB), .2);
                        }
                        dtlapusta = min(dtlapusta, dTETAmax * lrsf_temp / V);
                        
                        // Count old yelding nodes
                        bool ynn = false;
                        if (YNY0(i, j) > 0) {
                            ynn         = true;
                            DSY(i, j)   = SIIB(i, j) - syield;
                            ddd         += pow(DSY(i, j), 2);
                            ynpl++;
                        }

                        // Update viscosity
                        const double A = syield / siiel;
                        if (A < 1) {
                            // New viscosity for the basic node
                            etapl = dt * GGG(i, j) * A / (1 - A);
                            if (etapl < eta0_temp) {
                                // Update plastic nodes
                                ETA5(i, j)      = pow(etapl, 1 - etawt) * pow(ETA(i, j), etawt);
                                YNY(i, j)       = 1;
                                // Count yelding nodes
                                if (!ynn) {
                                    DSY(i, j)   = SIIB(i, j) - syield;
                                    ddd         += pow(DSY(i, j), 2);
                                    ynpl++;
                                }
                            } else {
                                ETA5(i, j) = eta0_temp;
                            }
                        } else {
                            ETA5(i, j) = eta0_temp;
                        }
                    }
                }
            }
        }
}