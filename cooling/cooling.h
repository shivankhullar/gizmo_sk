#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif

/*!
 * This file contains the definitions for the cooling.c routines
 */

/*!
 * This file was originally part of the GADGET3 code developed by Volker Springel.
 * The code has been modified by Phil Hopkins and Mike Grudic for GIZMO. Essentially everything has been re-written at this point.
 */

double ThermalProperties(double u, double rho, int target, double *mu_guess, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, struct particle_data *pp, struct gas_cell_data *cell);
double INLINE_FUNC return_uvb_shieldfac(int target, double gamma_12, double nHcgs, double logT, struct particle_data *pp, struct gas_cell_data *cell);
double INLINE_FUNC return_local_gammamultiplier(int target, struct particle_data *pp, struct gas_cell_data *cell);
double evaluate_Compton_heating_cooling_rate(int target, double T, double nHcgs, double n_elec, double shielding_factor_for_exgalbg, struct particle_data *pp, struct gas_cell_data *cell);
double get_background_radiation_temperature_for_emission_corrections(int target, struct particle_data *pp, struct gas_cell_data *cell);
void   InitCool(void);
#ifndef CHIMES 
void   InitCoolMemory(void);
void   IonizeParams(void);
void   IonizeParamsFunction(void);
void   IonizeParamsTable(void);
//double INLINE_FUNC LogTemp(double u, double ne);
void   MakeCoolingTable(void);
void   ReadIonizeParams(char *fname);
void   SetZeroIonization(void);
#endif 
void   TestCool(void);

#ifndef CHIMES 
double find_abundances_and_rates(double logT, double rho, int target, double shieldfac, int return_cooling_mode,
                                 double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess,
                                 double *LambdaExc_return, double *LambdaIon_return, double *LambdaRec_return, double *LambdaFF_return, struct particle_data *pp, struct gas_cell_data *cell);
double convert_u_to_temp(double u, double rho, int target, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess, struct particle_data *pp, struct gas_cell_data *cell);
double CoolingRate(double logT, double rho, double n_elec_guess, double *n_elec_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
double INLINE_FUNC CoolingRateFromU(double u, double rho, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
#endif
double DoCooling(double u_old, double rho, double dt, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
#ifndef CHIMES
double GetCoolingTime(double u_old, double rho, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
double DoInstabilityCooling(double m_old, double u, double rho, double dt, double fac, double ne_guess, double *ne_eval, int target, struct particle_data *pp, struct gas_cell_data *cell);
#endif

#ifdef COOL_GRACKLE
void InitGrackle(void);
double CallGrackle(double u_old, double rho, double dt, double ne_guess, int target, int mode);
#endif

