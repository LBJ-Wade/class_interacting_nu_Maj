/** @file background.c Documented background module
 *
 * * Julien Lesgourgues, 17.04.2011
 * * routines related to ncdm written by T. Tram in 2011
 *
 * Deals with the cosmological background evolution.
 * This module has two purposes:
 *
 * - at the beginning, to initialize the background, i.e. to integrate
 *    the background equations, and store all background quantities
 *    as a function of conformal time inside an interpolation table.
 *
 * - to provide routines which allow other modules to evaluate any
 *    background quantity for a given value of the conformal time (by
 *    interpolating within the interpolation table), or to find the
 *    correspondence between redshift and conformal time.
 *
 *
 * The overall logic in this module is the following:
 *
 * 1. most background parameters that we will call {A}
 * (e.g. rho_gamma, ..) can be expressed as simple analytical
 * functions of a few variables that we will call {B} (in simplest
 * models, of the scale factor 'a'; in extended cosmologies, of 'a'
 * plus e.g. (phi, phidot) for quintessence, or some temperature for
 * exotic particles, etc...).
 *
 * 2. in turn, quantities {B} can be found as a function of conformal
 * time by integrating the background equations.
 *
 * 3. some other quantities that we will call {C} (like e.g. the
 * sound horizon or proper time) also require an integration with
 * respect to time, that cannot be inferred analytically from
 * parameters {B}.
 *
 * So, we define the following routines:
 *
 * - background_functions() returns all background
 *    quantities {A} as a function of quantities {B}.
 *
 * - background_solve() integrates the quantities {B} and {C} with
 *    respect to conformal time; this integration requires many calls
 *    to background_functions().
 *
 * - the result is stored in the form of a big table in the background
 *    structure. There is one column for conformal time 'tau'; one or
 *    more for quantities {B}; then several columns for quantities {A}
 *    and {C}.
 *
 * Later in the code, if we know the variables {B} and need some
 * quantity {A}, the quickest and most precise way is to call directly
 * background_functions() (for instance, in simple models, if we want
 * H at a given value of the scale factor). If we know 'tau' and want
 * any other quantity, we can call background_at_tau(), which
 * interpolates in the table and returns all values. Finally it can be
 * useful to get 'tau' for a given redshift 'z': this can be done with
 * background_tau_of_z(). So if we are somewhere in the code, knowing
 * z and willing to get background quantities, we should call first
 * background_tau_of_z() and then background_at_tau().
 *
 *
 * In order to save time, background_at_tau() can be called in three
 * modes: short_info, normal_info, long_info (returning only essential
 * quantities, or useful quantities, or rarely useful
 * quantities). Each line in the interpolation table is a vector whose
 * first few elements correspond to the short_info format; a larger
 * fraction contribute to the normal format; and the full vector
 * corresponds to the long format. The guideline is that short_info
 * returns only geometric quantities like a, H, H'; normal format
 * returns quantities strictly needed at each step in the integration
 * of perturbations; long_info returns quantities needed only
 * occasionally.
 *
 * In summary, the following functions can be called from other modules:
 *
 * -# background_init() at the beginning
 * -# background_at_tau(), background_tau_of_z() at any later time
 * -# background_free() at the end, when no more calls to the previous functions are needed
 */

#include "background.h"
#include <stdbool.h>

/**
 * Background quantities at given conformal time tau.
 *
 * Evaluates all background quantities at a given value of
 * conformal time by reading the pre-computed table and interpolating.
 *
 * @param pba           Input: pointer to background structure (containing pre-computed table)
 * @param tau           Input: value of conformal time
 * @param return_format Input: format of output vector (short, normal, long)
 * @param intermode     Input: interpolation mode (normal or closeby)
 * @param last_index    Input/Output: index of the previous/current point in the interpolation array (input only for closeby mode, output for both)
 * @param pvecback      Output: vector (assumed to be already allocated)
 * @return the error status
 */

int background_at_tau(
                      struct background *pba,
                      double tau,
                      short return_format,
                      short intermode,
                      int * last_index,
                      double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                      ) {

  /** Summary: */

  /** - define local variables */

  /* size of output vector, controlled by input parameter return_format */
  int pvecback_size;

  /** - check that tau is in the pre-computed range */

  class_test(tau < pba->tau_table[0],
             pba->error_message,
             "out of range: tau=%e < tau_min=%e, you should decrease the precision parameter a_ini_over_a_today_default\n",tau,pba->tau_table[0]);

  class_test(tau > pba->tau_table[pba->bt_size-1],
             pba->error_message,
             "out of range: tau=%e > tau_max=%e\n",tau,pba->tau_table[pba->bt_size-1]);

  /** - deduce length of returned vector from format mode */

  if (return_format == pba->normal_info) {
    pvecback_size=pba->bg_size_normal;
  }
  else {
    if (return_format == pba->short_info) {
      pvecback_size=pba->bg_size_short;
    }
    else {
      pvecback_size=pba->bg_size;
    }
  }

  /** - interpolate from pre-computed table with array_interpolate()
      or array_interpolate_growing_closeby() (depending on
      interpolation mode) */

  if (intermode == pba->inter_normal) {
    class_call(array_interpolate_spline(
                                        pba->tau_table,
                                        pba->bt_size,
                                        pba->background_table,
                                        pba->d2background_dtau2_table,
                                        pba->bg_size,
                                        tau,
                                        last_index,
                                        pvecback,
                                        pvecback_size,
                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }
  if (intermode == pba->inter_closeby) {
    class_call(array_interpolate_spline_growing_closeby(
                                                        pba->tau_table,
                                                        pba->bt_size,
                                                        pba->background_table,
                                                        pba->d2background_dtau2_table,
                                                        pba->bg_size,
                                                        tau,
                                                        last_index,
                                                        pvecback,
                                                        pvecback_size,
                                                        pba->error_message),
               pba->error_message,
               pba->error_message);
  }

  return _SUCCESS_;
}

/**
 * Conformal time at given redshift.
 *
 * Returns tau(z) by interpolation from pre-computed table.
 *
 * @param pba Input: pointer to background structure
 * @param z   Input: redshift
 * @param tau Output: conformal time
 * @return the error status
 */

int background_tau_of_z(
                        struct background *pba,
                        double z,
                        double * tau
                        ) {

  /** Summary: */

  /** - define local variables */

  /* necessary for calling array_interpolate(), but never used */
  int last_index;

  /** - check that \f$ z \f$ is in the pre-computed range */
  class_test(z < pba->z_table[pba->bt_size-1],
             pba->error_message,
             "out of range: z=%e < z_min=%e\n",z,pba->z_table[pba->bt_size-1]);

  class_test(z > pba->z_table[0],
             pba->error_message,
             "out of range: a=%e > a_max=%e\n",z,pba->z_table[0]);

  /** - interpolate from pre-computed table with array_interpolate() */
  class_call(array_interpolate_spline(
                                      pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      pba->d2tau_dz2_table,
                                      1,
                                      z,
                                      &last_index,
                                      tau,
                                      1,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;
}

/**
 * Background quantities at given \f$ a \f$.
 *
 * Function evaluating all background quantities which can be computed
 * analytically as a function of {B} parameters such as the scale factor 'a'
 * (see discussion at the beginning of this file). In extended
 * cosmological models, the pvecback_B vector contains other input parameters than
 * just 'a', e.g. (phi, phidot) for quintessence, some temperature of
 * exotic relics, etc...
 *
 * @param pba           Input: pointer to background structure
 * @param pvecback_B    Input: vector containing all {B} type quantities (scale factor, ...)
 * @param return_format Input: format of output vector
 * @param pvecback      Output: vector of background quantities (assumed to be already allocated)
 * @return the error status
 */

int background_functions(
                         struct background *pba,
                         double * pvecback_B, /* Vector containing all {B} quantities. */
                         short return_format,
                         double * pvecback /* vector with argument pvecback[index_bg] (must be already allocated with a size compatible with return_format) */
                         ) {

  /** Summary: */

  /** - define local variables */

  /* total density */
  double rho_tot;
  /* critical density */
  double rho_crit;
  /* total pressure */
  double p_tot;
  /* total relativistic density */
  double rho_r;
  /* total non-relativistic density */
  double rho_m;
  /* scale factor relative to scale factor today */
  double a_rel;
  /* background ncdm quantities */
  double rho_ncdm,p_ncdm,pseudo_p_ncdm, T_ncdm, mu_ncdm, M, qmax;
  /* index for n_ncdm species */
  int n_ncdm, index_q;
  /* fluid's time-dependent equation of state parameter */
  double w_fld, dw_over_da, integral_fld;
  /* scale factor */
  double a;
  /* scalar field quantities */
  double phi, phi_prime;
  /* Since we only know a_prime_over_a after we have rho_tot,
     it is not possible to simply sum up p_tot_prime directly.
     Instead we sum up dp_dloga = p_prime/a_prime_over_a. The formula is
     p_prime = a_prime_over_a * dp_dloga = a_prime_over_a * Sum [ (w_prime/a_prime_over_a -3(1+w)w)rho].
     Note: The scalar field contribution must be added in the end, as an exception!*/
  double dp_dloga;

  /** - initialize local variables */
  a = pvecback_B[pba->index_bi_a];


  rho_tot = 0.;
  p_tot = 0.;
  dp_dloga = 0.;
  rho_r=0.;
  rho_m=0.;
  a_rel = a / pba->a_today;

  class_test(a_rel <= 0.,
             pba->error_message,
             "a = %e instead of strictly positive",a_rel);

  /** - pass value of \f$ a\f$ to output */
  pvecback[pba->index_bg_a] = a;

  /** - compute each component's density and pressure */

  /* photons */
  pvecback[pba->index_bg_rho_g] = pba->Omega0_g * pow(pba->H0,2) / pow(a_rel,4);
  rho_tot += pvecback[pba->index_bg_rho_g];
  p_tot += (1./3.) * pvecback[pba->index_bg_rho_g];
  dp_dloga += -(4./3.) * pvecback[pba->index_bg_rho_g];
  rho_r += pvecback[pba->index_bg_rho_g];

  /* baryons */
  pvecback[pba->index_bg_rho_b] = pba->Omega0_b * pow(pba->H0,2) / pow(a_rel,3);
  rho_tot += pvecback[pba->index_bg_rho_b];
  p_tot += 0;
  rho_m += pvecback[pba->index_bg_rho_b];

  /* cdm */
  if (pba->has_cdm == _TRUE_) {
    pvecback[pba->index_bg_rho_cdm] = pba->Omega0_cdm * pow(pba->H0,2) / pow(a_rel,3);
    rho_tot += pvecback[pba->index_bg_rho_cdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_cdm];
  }

  /* dcdm */
  if (pba->has_dcdm == _TRUE_) {
    /* Pass value of rho_dcdm to output */
    pvecback[pba->index_bg_rho_dcdm] = pvecback_B[pba->index_bi_rho_dcdm];
    rho_tot += pvecback[pba->index_bg_rho_dcdm];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_dcdm];
  }

  /* dr */
  if (pba->has_dr == _TRUE_) {
    /* Pass value of rho_dr to output */
    pvecback[pba->index_bg_rho_dr] = pvecback_B[pba->index_bi_rho_dr];
    rho_tot += pvecback[pba->index_bg_rho_dr];
    p_tot += (1./3.)*pvecback[pba->index_bg_rho_dr];
    dp_dloga += -(4./3.) * pvecback[pba->index_bg_rho_dr];
    rho_r += pvecback[pba->index_bg_rho_dr];
  }

  /* Scalar field */
  if (pba->has_scf == _TRUE_) {
    phi = pvecback_B[pba->index_bi_phi_scf];
    phi_prime = pvecback_B[pba->index_bi_phi_prime_scf];
    pvecback[pba->index_bg_phi_scf] = phi; // value of the scalar field phi
    pvecback[pba->index_bg_phi_prime_scf] = phi_prime; // value of the scalar field phi derivative wrt conformal time
    pvecback[pba->index_bg_V_scf] = V_scf(pba,phi); //V_scf(pba,phi); //write here potential as function of phi
    pvecback[pba->index_bg_dV_scf] = dV_scf(pba,phi); // dV_scf(pba,phi); //potential' as function of phi
    pvecback[pba->index_bg_ddV_scf] = ddV_scf(pba,phi); // ddV_scf(pba,phi); //potential'' as function of phi
    pvecback[pba->index_bg_rho_scf] = (phi_prime*phi_prime/(2*a*a) + V_scf(pba,phi))/3.; // energy of the scalar field. The field units are set automatically by setting the initial conditions
    pvecback[pba->index_bg_p_scf] =(phi_prime*phi_prime/(2*a*a) - V_scf(pba,phi))/3.; // pressure of the scalar field
    rho_tot += pvecback[pba->index_bg_rho_scf];
    p_tot += pvecback[pba->index_bg_p_scf];
    dp_dloga += 0.0; /** <-- This depends on a_prime_over_a, so we cannot add it now! */
    //divide relativistic & nonrelativistic (not very meaningful for oscillatory models)
    rho_r += 3.*pvecback[pba->index_bg_p_scf]; //field pressure contributes radiation
    rho_m += pvecback[pba->index_bg_rho_scf] - 3.* pvecback[pba->index_bg_p_scf]; //the rest contributes matter
    //printf(" a= %e, Omega_scf = %f, \n ",a_rel, pvecback[pba->index_bg_rho_scf]/rho_tot );
  }

  /* ncdm */
  if (pba->has_ncdm == _TRUE_) {

    /* Loop over species: */
    for(n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++){

      /* function returning background ncdm[n_ncdm] quantities (only
         those for which non-NULL pointers are passed) */

    if(pba->ncdm_background_distribution[n_ncdm] == _fermi_dirac_v2_ || pba->ncdm_background_distribution[n_ncdm] == _majoron_){
      class_call(interpolate_background_ncdm_distribution(pba,n_ncdm,pba->q_ncdm_bg[n_ncdm],pba->q_size_ncdm_bg[n_ncdm],1./a_rel-1.,pba->f_ncdm_bg[n_ncdm]),
      pba->error_message,
      pba->error_message);
      class_call(interpolate_T_and_mu_at_z(pba,n_ncdm,1./a_rel-1.,&T_ncdm,&mu_ncdm),
      pba->error_message,
      pba->error_message);
      // pvecback[pba->index_bg_T_ncdm1+n_ncdm] = T_ncdm*_eV_/_k_B_/pba->T_cmb;//Tncdm / Tcmb
      pvecback[pba->index_bg_T_ncdm1+n_ncdm] = T_ncdm;//Tncdm [eV]
      // pvecback[pba->index_bg_Mu_ncdm1+n_ncdm] = mu_ncdm*_eV_/_k_B_/pba->T_cmb;//mu_ncdm / Tcmb
      pvecback[pba->index_bg_Mu_ncdm1+n_ncdm] = mu_ncdm;//mu_ncdm [eV]
      // printf("n_ncdm %d mu_ncdm %e\n",n_ncdm,mu_ncdm);
      M = pba->m_ncdm_in_eV[n_ncdm];
      class_call(get_q_max(pba,n_ncdm,a_rel,M,&qmax),
      pba->error_message,
      pba->error_message);
    }else{
      for(index_q = 0; index_q < pba->q_size_ncdm_bg[n_ncdm]; index_q++){
      pba->f_ncdm_bg[n_ncdm][index_q] = 1;//f_ncdm is already included in w_ncdm_bg in the case of standard neutrinos.
      }
      qmax=1;
      M = pba->M_ncdm[n_ncdm];
      pvecback[pba->index_bg_T_ncdm1+n_ncdm] = pba->T_ncdm[n_ncdm]/a_rel*pba->T_cmb*_k_B_/_eV_; //converts to eV
      pvecback[pba->index_bg_Mu_ncdm1+n_ncdm] = pba->T_ncdm[n_ncdm]*pba->ksi_ncdm[n_ncdm]*pba->T_cmb*_k_B_/_eV_; //converts to eV
    }



      class_call(background_ncdm_momenta(
                                         pba->q_ncdm_bg[n_ncdm],
                                         pba->w_ncdm_bg[n_ncdm],
                                         pba->f_ncdm_bg[n_ncdm],
                                         pba->q_size_ncdm_bg[n_ncdm],
                                         M,
                                         qmax,
                                         pba->factor_ncdm[n_ncdm],
                                         1./a_rel-1.,
                                         n_ncdm,
                                         NULL,
                                         &rho_ncdm,
                                         &p_ncdm,
                                         NULL,
                                         &pseudo_p_ncdm),
                 pba->error_message,
                 pba->error_message);

      pvecback[pba->index_bg_rho_ncdm1+n_ncdm] = rho_ncdm;
      rho_tot += rho_ncdm;
      pvecback[pba->index_bg_p_ncdm1+n_ncdm] = p_ncdm;
      p_tot += p_ncdm;
      pvecback[pba->index_bg_pseudo_p_ncdm1+n_ncdm] = pseudo_p_ncdm;
      /** See e.g. Eq. A6 in 1811.00904. */
      dp_dloga += (pseudo_p_ncdm - 5*p_ncdm);

      /* (3 p_ncdm1) is the "relativistic" contribution to rho_ncdm1 */
      rho_r += 3.* p_ncdm;

      /* (rho_ncdm1 - 3 p_ncdm1) is the "non-relativistic" contribution
         to rho_ncdm1 */
      rho_m += rho_ncdm - 3.* p_ncdm;
    }
  }

  /* Lambda */
  if (pba->has_lambda == _TRUE_) {
    pvecback[pba->index_bg_rho_lambda] = pba->Omega0_lambda * pow(pba->H0,2);
    rho_tot += pvecback[pba->index_bg_rho_lambda];
    p_tot -= pvecback[pba->index_bg_rho_lambda];
  }

  /* fluid with w(a) and constant cs2 */
  if (pba->has_fld == _TRUE_) {

    /* get rho_fld from vector of integrated variables */
    pvecback[pba->index_bg_rho_fld] = pvecback_B[pba->index_bi_rho_fld];

    /* get w_fld from dedicated function */
    class_call(background_w_fld(pba,a,&w_fld,&dw_over_da,&integral_fld), pba->error_message, pba->error_message);
    pvecback[pba->index_bg_w_fld] = w_fld;

    // Obsolete: at the beginning, we had here the analytic integral solution corresponding to the case w=w0+w1(1-a/a0):
    // pvecback[pba->index_bg_rho_fld] = pba->Omega0_fld * pow(pba->H0,2) / pow(a_rel,3.*(1.+pba->w0_fld+pba->wa_fld)) * exp(3.*pba->wa_fld*(a_rel-1.));
    // But now everthing is integrated numerically for a given w_fld(a) defined in the function background_w_fld.

    rho_tot += pvecback[pba->index_bg_rho_fld];
    p_tot += w_fld * pvecback[pba->index_bg_rho_fld];
    dp_dloga += (a*dw_over_da-3*(1+w_fld)*w_fld)*pvecback[pba->index_bg_rho_fld];
  }

  /* relativistic neutrinos (and all relativistic relics) */
  if (pba->has_ur == _TRUE_) {
    pvecback[pba->index_bg_rho_ur] = pba->Omega0_ur * pow(pba->H0,2) / pow(a_rel,4);
    rho_tot += pvecback[pba->index_bg_rho_ur];
    p_tot += (1./3.) * pvecback[pba->index_bg_rho_ur];
    dp_dloga += -(4./3.) * pvecback[pba->index_bg_rho_ur];
    rho_r += pvecback[pba->index_bg_rho_ur];
  }

  /* interacting dark matter */
  if (pba->has_idm_dr == _TRUE_) {
    pvecback[pba->index_bg_rho_idm_dr] = pba->Omega0_idm_dr * pow(pba->H0,2) / pow(a_rel,3);
    rho_tot += pvecback[pba->index_bg_rho_idm_dr];
    p_tot += 0.;
    rho_m += pvecback[pba->index_bg_rho_idm_dr];
  }

  /* interacting dark radiation */
  if (pba->has_idr == _TRUE_) {
    pvecback[pba->index_bg_rho_idr] = pba->Omega0_idr * pow(pba->H0,2) / pow(a_rel,4);
    rho_tot += pvecback[pba->index_bg_rho_idr];
    p_tot += (1./3.) * pvecback[pba->index_bg_rho_idr];
    rho_r += pvecback[pba->index_bg_rho_idr];
  }

  /** - compute expansion rate H from Friedmann equation: this is the
      only place where the Friedmann equation is assumed. Remember
      that densities are all expressed in units of \f$ [3c^2/8\pi G] \f$, ie
      \f$ \rho_{class} = [8 \pi G \rho_{physical} / 3 c^2]\f$ */
  pvecback[pba->index_bg_H] = sqrt(rho_tot-pba->K/a/a);

  /** - compute derivative of H with respect to conformal time */
  pvecback[pba->index_bg_H_prime] = - (3./2.) * (rho_tot + p_tot) * a + pba->K/a;

  /* Total energy density*/
  pvecback[pba->index_bg_rho_tot] = rho_tot;

  /* Total pressure */
  pvecback[pba->index_bg_p_tot] = p_tot;

  /* Derivative of total pressure w.r.t. conformal time */
  pvecback[pba->index_bg_p_tot_prime] = a*pvecback[pba->index_bg_H]*dp_dloga;
  if (pba->has_scf == _TRUE_){
    /** The contribution of scf was not added to dp_dloga, add p_scf_prime here: */
    pvecback[pba->index_bg_p_prime_scf] = pvecback[pba->index_bg_phi_prime_scf]*
      (-pvecback[pba->index_bg_phi_prime_scf]*pvecback[pba->index_bg_H]/a-2./3.*pvecback[pba->index_bg_dV_scf]);
    pvecback[pba->index_bg_p_tot_prime] += pvecback[pba->index_bg_p_prime_scf];
  }

  /** - compute critical density */
  rho_crit = rho_tot-pba->K/a/a;
  class_test(rho_crit <= 0.,
             pba->error_message,
             "rho_crit = %e instead of strictly positive",rho_crit);

  /** - compute relativistic density to total density ratio */
  pvecback[pba->index_bg_Omega_r] = rho_r / rho_crit;

  /** - compute other quantities in the exhaustive, redundant format */
  if (return_format == pba->long_info) {

    /** - store critical density */
    pvecback[pba->index_bg_rho_crit] = rho_crit;

    /** - compute Omega_m */
    pvecback[pba->index_bg_Omega_m] = rho_m / rho_crit;

    /* one can put other variables here */
    /*  */
    /*  */

  }

  return _SUCCESS_;

}

/**
 * Single place where the fluid equation of state is
 * defined. Parameters of the function are passed through the
 * background structure. Generalisation to arbitrary functions should
 * be simple.
 *
 * @param pba            Input: pointer to background structure
 * @param a              Input: current value of scale factor
 * @param w_fld          Output: equation of state parameter w_fld(a)
 * @param dw_over_da_fld Output: function dw_fld/da
 * @param integral_fld   Output: function \f$ \int_{a}^{a_0} da 3(1+w_{fld})/a \f$
 * @return the error status
 */

int background_w_fld(
                     struct background * pba,
                     double a,
                     double * w_fld,
                     double * dw_over_da_fld,
                     double * integral_fld) {

  double Omega_ede = 0.;
  double dOmega_ede_over_da = 0.;
  double d2Omega_ede_over_da2 = 0.;
  double a_eq, Omega_r, Omega_m;

  /** - first, define the function w(a) */
  switch (pba->fluid_equation_of_state) {
  case CLP:
    *w_fld = pba->w0_fld + pba->wa_fld * (1. - a / pba->a_today);
    break;
  case EDE:
    // Omega_ede(a) taken from eq. (10) in 1706.00730
    Omega_ede = (pba->Omega0_fld - pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld)))
      /(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld))
      + pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld));

    // d Omega_ede / d a taken analytically from the above
    dOmega_ede_over_da = - pba->Omega_EDE* 3.*pba->w0_fld*pow(a,-3.*pba->w0_fld-1.)/(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld))
      - (pba->Omega0_fld - pba->Omega_EDE*(1.-pow(a,-3.*pba->w0_fld)))*(1.-pba->Omega0_fld)*3.*pba->w0_fld*pow(a,3.*pba->w0_fld-1.)/pow(pba->Omega0_fld+(1.-pba->Omega0_fld)*pow(a,3.*pba->w0_fld),2)
      + pba->Omega_EDE*3.*pba->w0_fld*pow(a,-3.*pba->w0_fld-1.);

    // find a_equality (needed because EDE tracks first radiation, then matter)
    Omega_r = pba->Omega0_g * (1. + 3.046 * 7./8.*pow(4./11.,4./3.)); // assumes LambdaCDM + eventually massive neutrinos so light that they are relativistic at equality; needs to be generalised later on.
    Omega_m = pba->Omega0_b;
    if (pba->has_cdm == _TRUE_) Omega_m += pba->Omega0_cdm;
    if (pba->has_idm_dr == _TRUE_) Omega_m += pba->Omega0_idm_dr;
    if (pba->has_dcdm == _TRUE_)
      class_stop(pba->error_message,"Early Dark Energy not compatible with decaying Dark Matter because we omitted to code the calculation of a_eq in that case, but it would not be difficult to add it if necessary, should be a matter of 5 minutes");
    a_eq = Omega_r/Omega_m; // assumes a flat universe with a=1 today

    // w_ede(a) taken from eq. (11) in 1706.00730
    *w_fld = - dOmega_ede_over_da*a/Omega_ede/3./(1.-Omega_ede)+a_eq/3./(a+a_eq);
    break;
  }


  /** - then, give the corresponding analytic derivative dw/da (used
      by perturbation equations; we could compute it numerically,
      but with a loss of precision; as long as there is a simple
      analytic expression of the derivative of the previous
      function, let's use it! */
  switch (pba->fluid_equation_of_state) {
  case CLP:
    *dw_over_da_fld = - pba->wa_fld / pba->a_today;
    break;
  case EDE:
    d2Omega_ede_over_da2 = 0.;
    *dw_over_da_fld = - d2Omega_ede_over_da2*a/3./(1.-Omega_ede)/Omega_ede
      - dOmega_ede_over_da/3./(1.-Omega_ede)/Omega_ede
      + dOmega_ede_over_da*dOmega_ede_over_da*a/3./(1.-Omega_ede)/(1.-Omega_ede)/Omega_ede
      + a_eq/3./(a+a_eq)/(a+a_eq);
    break;
  }

  /** - finally, give the analytic solution of the following integral:
        \f$ \int_{a}^{a0} da 3(1+w_{fld})/a \f$. This is used in only
        one place, in the initial conditions for the background, and
        with a=a_ini. If your w(a) does not lead to a simple analytic
        solution of this integral, no worry: instead of writing
        something here, the best would then be to leave it equal to
        zero, and then in background_initial_conditions() you should
        implement a numerical calculation of this integral only for
        a=a_ini, using for instance Romberg integration. It should be
        fast, simple, and accurate enough. */
  switch (pba->fluid_equation_of_state) {
  case CLP:
    *integral_fld = 3.*((1.+pba->w0_fld+pba->wa_fld)*log(pba->a_today/a) + pba->wa_fld*(a/pba->a_today-1.));
    break;
  case EDE:
    class_stop(pba->error_message,"EDE implementation not finished: to finish it, read the comments in background.c just before this line\n");
    break;
  }

  /** note: of course you can generalise these formulas to anything,
      defining new parameters pba->w..._fld. Just remember that so
      far, HyRec explicitely assumes that w(a)= w0 + wa (1-a/a0); but
      Recfast does not assume anything */

  return _SUCCESS_;
}

/**
 * Initialize the background structure, and in particular the
 * background interpolation table.
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input/Output: pointer to initialized background structure
 * @return the error status
 */

int background_init(
                    struct precision * ppr,
                    struct background * pba
                    ) {

  /** Summary: */

  /** - define local variables */
  int n_ncdm, index_q;
  double rho_ncdm_rel,rho_nu_rel, M, T_ncdm, mu_ncdm, qmax;
  double Neff, N_dark;
  double w_fld, dw_over_da, integral_fld;
  int filenum=0;

  /** - in verbose mode, provide some information */
  if (pba->background_verbose > 0) {
    printf("Running CLASS version %s\n",_VERSION_);
    printf("Computing background\n");

    /* below we want to inform the user about ncdm species and/or the total N_eff */
    if ((pba->N_ncdm > 0) || (pba->Omega0_idr != 0.))  {

      /* contribution of ultra-relativistic species _ur to N_eff */
      Neff = pba->Omega0_ur/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;

      /* contribution of ncdm species to N_eff*/
      if (pba->N_ncdm > 0){
        /* loop over ncdm species */
        for (n_ncdm=0;n_ncdm<pba->N_ncdm; n_ncdm++) {

          /* inform if p-s-d read in files */
          if (pba->got_files[n_ncdm] == _TRUE_) {
            printf(" -> ncdm species i=%d read from file %s\n",n_ncdm+1,pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
            filenum++;
          }

        if(pba->ncdm_background_distribution[n_ncdm] == _fermi_dirac_v2_ || pba->ncdm_background_distribution[n_ncdm] == _majoron_){
            // printf("1./(ppr->a_ini_over_a_today_default * pba->a_today)-1. %e\n", 1./(ppr->a_ini_over_a_today_default * pba->a_today)-1.);
            class_call(interpolate_background_ncdm_distribution(pba,n_ncdm,pba->q_ncdm_bg[n_ncdm],pba->q_size_ncdm_bg[n_ncdm],1./(ppr->a_ini_over_a_today_default * pba->a_today)-1.,pba->f_ncdm_bg[n_ncdm]),
            // class_call(interpolate_background_ncdm_distribution(pba,n_ncdm,0.),
            pba->error_message,
            pba->error_message);

            M = pba->m_ncdm_in_eV[n_ncdm];
            class_call(get_q_max(pba,n_ncdm,ppr->a_ini_over_a_today_default * pba->a_today,M,&qmax),
            // class_call(get_q_max(pba,n_ncdm,1.,M,&qmax),
            pba->error_message,
            pba->error_message);

            // printf("right after bug\n");
          }else{
            for(index_q = 0; index_q < pba->q_size_ncdm_bg[n_ncdm]; index_q++){
            pba->f_ncdm_bg[n_ncdm][index_q] = 1;//f_ncdm is already included in w_ncdm_bg in the case of standard neutrinos.
            }
            qmax = 1; //qmax alrady taken into account when defning integral weights.
            M = pba->M_ncdm[n_ncdm];
          }


          /* call this function to get rho_ncdm */
          background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                  pba->w_ncdm_bg[n_ncdm],
                                  pba->f_ncdm_bg[n_ncdm],
                                  pba->q_size_ncdm_bg[n_ncdm],
                                  // pba->M_ncdm[n_ncdm],
                                  M,
                                  qmax,
                                  pba->factor_ncdm[n_ncdm],
                                  1./(ppr->a_ini_over_a_today_default * pba->a_today)-1.,
                                  n_ncdm,
                                  NULL,
                                  &rho_ncdm_rel,
                                  NULL,
                                  NULL,
                                  NULL);

          /* inform user of the contribution of each species to
             radiation density (in relativistic limit): should be
             between 1.01 and 1.02 for each active neutrino species;
             evaluated as rho_ncdm/rho_nu_rel where rho_nu_rel is the
             density of one neutrino in the instantaneous decoupling
             limit, i.e. assuming T_nu=(4/11)^1/3 T_gamma (this comes
             from the definition of N_eff) */
          rho_nu_rel = 56.0/45.0*pow(_PI_,6)*pow(4.0/11.0,4.0/3.0)*_G_/pow(_h_P_,3)/pow(_c_,7)*
            pow(_Mpc_over_m_,2)*pow(pba->T_cmb*_k_B_,4)/(pow(ppr->a_ini_over_a_today_default * pba->a_today,4));

          printf(" -> ncdm species i=%d sampled with %d (resp. %d) points for purpose of background (resp. perturbation) integration. In the relativistic limit:rho = %e & rho_nu %e, it gives Delta N_eff = %g\n",
                 n_ncdm+1,
                 pba->q_size_ncdm_bg[n_ncdm],
                 pba->q_size_ncdm[n_ncdm],
                 rho_ncdm_rel,
                 rho_nu_rel,
                 rho_ncdm_rel/rho_nu_rel);

          Neff += rho_ncdm_rel/rho_nu_rel;
        }
      }

      /* contribution of interacting dark radiation _idr to N_eff */
      if (pba->Omega0_idr != 0.) {
        N_dark = pba->Omega0_idr/7.*8./pow(4./11.,4./3.)/pba->Omega0_g;
        Neff += N_dark;
        printf(" -> dark radiation Delta Neff %e\n",N_dark);
      }

      printf(" -> total N_eff = %g (sumed over ultra-relativistic species, ncdm and dark radiation)\n",Neff);

    }
  }

  /** - if shooting failed during input, catch the error here */
  class_test_except(pba->shooting_failed == _TRUE_,
                    pba->error_message,
                    background_free_input(pba),
                    "Shooting failed, try optimising input_get_guess(). Error message:\n\n%s",
                    pba->shooting_error);

  /** - assign values to all indices in vectors of background quantities with background_indices()*/
  class_call(background_indices(pba),
             pba->error_message,
             pba->error_message);

  /* fluid equation of state */
  if (pba->has_fld == _TRUE_) {

    class_call(background_w_fld(pba,0.,&w_fld,&dw_over_da,&integral_fld), pba->error_message, pba->error_message);

    class_test(w_fld >= 1./3.,
               pba->error_message,
               "Your choice for w(a--->0)=%g is suspicious, since it is bigger than -1/3 there cannot be radiation domination at early times\n",
               w_fld);
  }

  /* in verbose mode, inform the user about the value of the ncdm
     masses in eV and about the ratio [m/omega_ncdm] in eV (the usual
     93 point something)*/
  if ((pba->background_verbose > 0) && (pba->has_ncdm == _TRUE_)) {
    for (n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
      printf(" -> non-cold dark matter species with i=%d has m_i = %e eV (so m_i / omega_i =%e eV)\n",
             n_ncdm+1,
             pba->m_ncdm_in_eV[n_ncdm],
             pba->m_ncdm_in_eV[n_ncdm]*pba->deg_ncdm[n_ncdm]/pba->Omega0_ncdm[n_ncdm]/pba->h/pba->h);
    }
  }

  /* check other quantities which would lead to segmentation fault if zero */
  class_test(pba->a_today <= 0,
             pba->error_message,
             "input a_today = %e instead of strictly positive",pba->a_today);

  class_test(_Gyr_over_Mpc_ <= 0,
             pba->error_message,
             "_Gyr_over_Mpc = %e instead of strictly positive",_Gyr_over_Mpc_);

  /** - this function integrates the background over time, allocates
      and fills the background table */
  class_call(background_solve(ppr,pba),
             pba->error_message,
             pba->error_message);

  /** - this function finds and stores a few derived parameters at radiation-matter equality */
  class_call(background_find_equality(ppr,pba),
             pba->error_message,
             pba->error_message);

  class_call(background_output_budget(pba),
             pba->error_message,
             pba->error_message);

  return _SUCCESS_;

}

/**
 * Free all memory space allocated by background_init().
 *
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */

int background_free(
                    struct background *pba
                    ) {

  class_call(background_free_noinput(pba),
              pba->error_message,
              pba->error_message);

  class_call(background_free_input(pba),
              pba->error_message,
              pba->error_message);

  return _SUCCESS_;
}

/**
 * Free only the memory space NOT allocated through input_read_parameters()
 *
 * @param pba Input: pointer to background structure (to be freed)
 * @return the error status
 */

int background_free_noinput(
                            struct background *pba
                            ) {

  free(pba->tau_table);
  free(pba->z_table);
  free(pba->d2tau_dz2_table);
  free(pba->background_table);
  free(pba->d2background_dtau2_table);

  return _SUCCESS_;
}
/**
 * Free pointers inside background structure which were
 * allocated in input_read_parameters()
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */

int background_free_input(
                          struct background *pba
                          ) {

  int k;

  if (pba->Omega0_ncdm_tot != 0.){
    for(k=0; k<pba->N_ncdm; k++){
      free(pba->q_ncdm[k]);
      free(pba->w_ncdm[k]);
      free(pba->q_ncdm_bg[k]);
      free(pba->w_ncdm_bg[k]);
      free(pba->dlnf0_dlnq_ncdm[k]);
      free(pba->f_ncdm_bg[k]);
      free(pba->f_ncdm[k]);
    }
    free(pba->ncdm_quadrature_strategy);
    free(pba->ncdm_input_q_size_bg);
    free(pba->ncdm_input_q_size);
    free(pba->ncdm_qmax);
    free(pba->q_ncdm);
    free(pba->w_ncdm);
    free(pba->q_ncdm_bg);
    free(pba->w_ncdm_bg);
    free(pba->dlnf0_dlnq_ncdm);
    free(pba->f_ncdm_bg);
    free(pba->f_ncdm);
    free(pba->q_size_ncdm);
    free(pba->q_size_ncdm_bg);
    free(pba->M_ncdm);
    free(pba->T_ncdm);
    free(pba->ksi_ncdm);
    free(pba->deg_ncdm);
    free(pba->Omega0_ncdm);
    free(pba->m_ncdm_in_eV);
    free(pba->z_nrel);
    free(pba->factor_ncdm);
    if(pba->got_files!=NULL)
      free(pba->got_files);
    if(pba->ncdm_psd_files!=NULL)
      free(pba->ncdm_psd_files);
    if(pba->ncdm_psd_parameters!=NULL)
      free(pba->ncdm_psd_parameters);
  }

  if (pba->Omega0_scf != 0.){
    if (pba->scf_parameters != NULL)
      free(pba->scf_parameters);
  }
  return _SUCCESS_;
}

/**
 * Assign value to each relevant index in vectors of background quantities.
 *
 * @param pba Input: pointer to background structure
 * @return the error status
 */

int background_indices(
                       struct background *pba
                       ) {

  /** Summary: */

  /** - define local variables */

  /* a running index for the vector of background quantities */
  int index_bg;
  /* a running index for the vector of background quantities to be integrated */
  int index_bi;

  /** - initialize all flags: which species are present? */

  pba->has_cdm = _FALSE_;
  pba->has_ncdm = _FALSE_;
  pba->has_dcdm = _FALSE_;
  pba->has_dr = _FALSE_;
  pba->has_scf = _FALSE_;
  pba->has_lambda = _FALSE_;
  pba->has_fld = _FALSE_;
  pba->has_ur = _FALSE_;
  pba->has_idr = _FALSE_;
  pba->has_idm_dr = _FALSE_;
  pba->has_curvature = _FALSE_;

  if (pba->Omega0_cdm != 0.)
    pba->has_cdm = _TRUE_;

  if (pba->Omega0_ncdm_tot != 0.)
    pba->has_ncdm = _TRUE_;

  if (pba->Omega0_dcdmdr != 0.){
    pba->has_dcdm = _TRUE_;
    if (pba->Gamma_dcdm != 0.)
      pba->has_dr = _TRUE_;
  }

  if (pba->Omega0_scf != 0.)
    pba->has_scf = _TRUE_;

  if (pba->Omega0_lambda != 0.)
    pba->has_lambda = _TRUE_;

  if (pba->Omega0_fld != 0.)
    pba->has_fld = _TRUE_;

  if (pba->Omega0_ur != 0.)
    pba->has_ur = _TRUE_;

  if (pba->Omega0_idr != 0.)
    pba->has_idr = _TRUE_;

  if (pba->Omega0_idm_dr != 0.)
    pba->has_idm_dr = _TRUE_;

  if (pba->sgnK != 0)
    pba->has_curvature = _TRUE_;

  /** - initialize all indices */

  index_bg=0;

  /* index for scale factor */
  class_define_index(pba->index_bg_a,_TRUE_,index_bg,1);

  /* - indices for H and its conformal-time-derivative */
  class_define_index(pba->index_bg_H,_TRUE_,index_bg,1);
  class_define_index(pba->index_bg_H_prime,_TRUE_,index_bg,1);

  /* - end of indices in the short vector of background values */
  pba->bg_size_short = index_bg;

  /* - index for rho_g (photon density) */
  class_define_index(pba->index_bg_rho_g,_TRUE_,index_bg,1);

  /* - index for rho_b (baryon density) */
  class_define_index(pba->index_bg_rho_b,_TRUE_,index_bg,1);

  /* - index for rho_cdm */
  class_define_index(pba->index_bg_rho_cdm,pba->has_cdm,index_bg,1);

  /* - indices for ncdm. We only define the indices for ncdm1
     (density, pressure, pseudo-pressure), the other ncdm indices
     are contiguous */
  class_define_index(pba->index_bg_rho_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_pseudo_p_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_T_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_dT_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_ddT_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_Mu_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_dMu_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);
  class_define_index(pba->index_bg_ddMu_ncdm1,pba->has_ncdm,index_bg,pba->N_ncdm);

  /* - index for dcdm */
  class_define_index(pba->index_bg_rho_dcdm,pba->has_dcdm,index_bg,1);

  /* - index for dr */
  class_define_index(pba->index_bg_rho_dr,pba->has_dr,index_bg,1);

  /* - indices for scalar field */
  class_define_index(pba->index_bg_phi_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_phi_prime_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_V_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_dV_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_ddV_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_rho_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_p_scf,pba->has_scf,index_bg,1);
  class_define_index(pba->index_bg_p_prime_scf,pba->has_scf,index_bg,1);

  /* - index for Lambda */
  class_define_index(pba->index_bg_rho_lambda,pba->has_lambda,index_bg,1);

  /* - index for fluid */
  class_define_index(pba->index_bg_rho_fld,pba->has_fld,index_bg,1);
  class_define_index(pba->index_bg_w_fld,pba->has_fld,index_bg,1);

  /* - index for ultra-relativistic neutrinos/species */
  class_define_index(pba->index_bg_rho_ur,pba->has_ur,index_bg,1);

  /* - index for total density */
  class_define_index(pba->index_bg_rho_tot,_TRUE_,index_bg,1);

  /* - index for total pressure */
  class_define_index(pba->index_bg_p_tot,_TRUE_,index_bg,1);

  /* - index for derivative of total pressure */
  class_define_index(pba->index_bg_p_tot_prime,_TRUE_,index_bg,1);

  /* - index for Omega_r (relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_r,_TRUE_,index_bg,1);

  /* - index interacting for dark radiation */
  class_define_index(pba->index_bg_rho_idr,pba->has_idr,index_bg,1);

  /* - index for interacting dark matter */
  class_define_index(pba->index_bg_rho_idm_dr,pba->has_idm_dr,index_bg,1);

  /* - put here additional ingredients that you want to appear in the
     normal vector */
  /*    */
  /*    */

  /* - end of indices in the normal vector of background values */
  pba->bg_size_normal = index_bg;

  /* - indices in the long version : */

  /* -> critical density */
  class_define_index(pba->index_bg_rho_crit,_TRUE_,index_bg,1);

  /* - index for Omega_m (non-relativistic density fraction) */
  class_define_index(pba->index_bg_Omega_m,_TRUE_,index_bg,1);

  /* -> conformal distance */
  class_define_index(pba->index_bg_conf_distance,_TRUE_,index_bg,1);

  /* -> angular diameter distance */
  class_define_index(pba->index_bg_ang_distance,_TRUE_,index_bg,1);

  /* -> luminosity distance */
  class_define_index(pba->index_bg_lum_distance,_TRUE_,index_bg,1);

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bg_time,_TRUE_,index_bg,1);

  /* -> conformal sound horizon */
  class_define_index(pba->index_bg_rs,_TRUE_,index_bg,1);

  /* -> density growth factor in dust universe */
  class_define_index(pba->index_bg_D,_TRUE_,index_bg,1);

  /* -> velocity growth factor in dust universe */
  class_define_index(pba->index_bg_f,_TRUE_,index_bg,1);

  /* -> put here additional quantities describing background */
  /*    */
  /*    */

  /* -> end of indices in the long vector of background values */
  pba->bg_size = index_bg;

  /* - now, indices in vector of variables to integrate.
     First {B} variables, then {C} variables. */

  index_bi=0;

  /* -> scale factor */
  class_define_index(pba->index_bi_a,_TRUE_,index_bi,1);

  /* -> energy density in DCDM */
  class_define_index(pba->index_bi_rho_dcdm,pba->has_dcdm,index_bi,1);

  /* -> energy density in DR */
  class_define_index(pba->index_bi_rho_dr,pba->has_dr,index_bi,1);

  /* -> energy density in fluid */
  class_define_index(pba->index_bi_rho_fld,pba->has_fld,index_bi,1);

  /* -> scalar field and its derivative wrt conformal time (Zuma) */
  class_define_index(pba->index_bi_phi_scf,pba->has_scf,index_bi,1);
  class_define_index(pba->index_bi_phi_prime_scf,pba->has_scf,index_bi,1);

  /* End of {B} variables, now continue with {C} variables */
  pba->bi_B_size = index_bi;

  /* -> proper time (for age of the Universe) */
  class_define_index(pba->index_bi_time,_TRUE_,index_bi,1);

  /* -> sound horizon */
  class_define_index(pba->index_bi_rs,_TRUE_,index_bi,1);

  /* -> Second order equation for growth factor */
  class_define_index(pba->index_bi_D,_TRUE_,index_bi,1);
  class_define_index(pba->index_bi_D_prime,_TRUE_,index_bi,1);

  /* -> index for conformal time in vector of variables to integrate */
  class_define_index(pba->index_bi_tau,_TRUE_,index_bi,1);

  /* -> end of indices in the vector of variables to integrate */
  pba->bi_size = index_bi;

  /* index_bi_tau must be the last index, because tau is part of this vector for the purpose of being stored, */
  /* but it is not a quantity to be integrated (since integration is over tau itself) */
  class_test(pba->index_bi_tau != index_bi-1,
             pba->error_message,
             "background integration requires index_bi_tau to be the last of all index_bi's");

  /* flags for calling the interpolation routine */

  pba->short_info=0;
  pba->normal_info=1;
  pba->long_info=2;

  pba->inter_normal=0;
  pba->inter_closeby=1;

  return _SUCCESS_;

}

/**
 * This is the routine where the distribution function f0(q) of each
 * ncdm species is specified (it is the only place to modify if you
 * need a partlar f0(q))
 *
 * @param pbadist Input:  structure containing all parameters defining f0(q)
 * @param q       Input:  momentum
 * @param f0      Output: phase-space distribution
 */

int background_ncdm_distribution(
                                 void * pbadist,
                                 double q,
                                 double * f0,
                                 double z
                                 ) {
  struct background * pba;
  struct background_parameters_for_distributions * pbadist_local;
  int n_ncdm,lastidx,last_index;
  double T_ncdm,mu_ncdm,T0_ncdm;
  double ksi,qmax;
  double qlast,dqlast,f0last,df0last;
  double *param;
  double eps;

  /* Variables corresponding to entries in param: */
  //double square_s12,square_s23,square_s13;
  //double mixing_matrix[3][3];
  //int i;

  /** - extract from the input structure pbadist all the relevant information */
  pbadist_local = pbadist;          /* restore actual format of pbadist */
  pba = pbadist_local->pba;         /* extract the background structure from it */
  param = pba->ncdm_psd_parameters; /* extract the optional parameter list from it */
  n_ncdm = pbadist_local->n_ncdm;   /* extract index of ncdm species under consideration */
  ksi = pba->ksi_ncdm[n_ncdm];      /* extract chemical potential */

  /** - shall we interpolate in file, or shall we use analytical formula below? */

//   SJW
// Note: redshifts (z_maj), Temp Majoron [eV] (T_maj), Temp Nu [eV] (T_nu), Chem Pot Maj [eV] (Mu_maj), Chem Pot Nu [eV] (Mu_nu)
// length of array: pba->len_maj


//      This is a quick test....
//
//      for (int ii=0; ii <= pba->len_maj; ii++){
//        printf("HERE: %d \t %e \t %e \t %e \t %e \t %e \n", ii, pba->z_maj[ii], pba->T_maj[ii], pba->T_nu[ii], pba->Mu_maj[ii], pba->Mu_nu[ii]);
//    }
//    exit(0);

  /** - a) deal first with the case of interpolating in files */
  if (pba->got_files[n_ncdm]==_TRUE_) {

    lastidx = pbadist_local->tablesize-1;
    if(q<pbadist_local->q[0]){
      //Handle q->0 case:
      *f0 = pbadist_local->f0[0];
    }
    else if(q>pbadist_local->q[lastidx]){
      //Handle q>qmax case (ensure continuous and derivable function with Boltzmann tail):
      qlast=pbadist_local->q[lastidx];
      f0last=pbadist_local->f0[lastidx];
      dqlast=qlast - pbadist_local->q[lastidx-1];
      df0last=f0last - pbadist_local->f0[lastidx-1];

      *f0 = f0last*exp(-(qlast-q)*df0last/f0last/dqlast);
    }
    else{
      //Do interpolation:
      class_call(array_interpolate_spline(
                                          pbadist_local->q,
                                          pbadist_local->tablesize,
                                          pbadist_local->f0,
                                          pbadist_local->d2f0,
                                          1,
                                          q,
                                          &pbadist_local->last_index,
                                          f0,
                                          1,
                                          pba->error_message),
                 pba->error_message,     pba->error_message);
    }
  }

  /** - b) deal now with case of reading analytical function */
  else{
    /**
       Next enter your analytic expression(s) for the p.s.d.'s. If
       you need different p.s.d.'s for different species, put each
       p.s.d inside a condition, like for instance: if (n_ncdm==2)
       {*f0=...}.  Remember that n_ncdm = 0 refers to the first
       species.
    */

    if(pba->ncdm_background_distribution[n_ncdm]==_fermi_dirac_){

      //VP: the standard neutrino distribution from CLASS


      /**************************************************/
      /*    FERMI-DIRAC INCLUDING CHEMICAL POTENTIALS   */
      /**************************************************/

      *f0 = 1.0/pow(2*_PI_,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));

      /**************************************************/

      /** This form is only appropriate for approximate studies, since in
          reality the chemical potentials are associated with flavor
          eigenstates, not mass eigenstates. It is easy to take this into
          account by introducing the mixing angles. In the later part
          (not read by the code) we illustrate how to do this. */

      if (_FALSE_) {

        /* We must use the list of extra parameters read in input, stored in the
           ncdm_psd_parameter list, extracted above from the structure
           and now called param[..] */

        /* check that this list has been read */
        class_test(param == NULL,
                   pba->error_message,
                   "Analytic expression wants to use 'ncdm_psd_parameters', but they have not been entered!");

        /* extract values from the list (in this example, mixing angles) */
        double square_s12=param[0];
        double square_s23=param[1];
        double square_s13=param[2];

        /* infer mixing matrix */
        double mixing_matrix[3][3];
        int i;

        mixing_matrix[0][0]=pow(fabs(sqrt((1-square_s12)*(1-square_s13))),2);
        mixing_matrix[0][1]=pow(fabs(sqrt(square_s12*(1-square_s13))),2);
        mixing_matrix[0][2]=fabs(square_s13);
        mixing_matrix[1][0]=pow(fabs(sqrt((1-square_s12)*square_s13*square_s23)+sqrt(square_s12*(1-square_s23))),2);
        mixing_matrix[1][1]=pow(fabs(sqrt(square_s12*square_s23*square_s13)-sqrt((1-square_s12)*(1-square_s23))),2);
        mixing_matrix[1][2]=pow(fabs(sqrt(square_s23*(1-square_s13))),2);
        mixing_matrix[2][0]=pow(fabs(sqrt(square_s12*square_s23)-sqrt((1-square_s12)*square_s13*(1-square_s23))),2);
        mixing_matrix[2][1]=pow(sqrt((1-square_s12)*square_s23)+sqrt(square_s12*square_s13*(1-square_s23)),2);
        mixing_matrix[2][2]=pow(fabs(sqrt((1-square_s13)*(1-square_s23))),2);

        /* loop over flavor eigenstates and compute psd of mass eigenstates */
        *f0=0.0;
        for(i=0;i<3;i++){

      	*f0 += mixing_matrix[i][n_ncdm]*1.0/pow(2*_PI_,3)*(1./(exp(q-pba->ksi_ncdm[i])+1.) +1./(exp(q+pba->ksi_ncdm[i])+1.));

        }
      } /* end of region not used, but shown as an example */
    }else if(pba->ncdm_background_distribution[n_ncdm]==_fermi_dirac_v2_ || pba->ncdm_background_distribution[n_ncdm]==_majoron_){
      //VP: here we define the background distribution
      interpolate_T_and_mu_at_z(pba,n_ncdm,z,&T_ncdm,&mu_ncdm);

      class_call(get_q_max(pba,n_ncdm,1./(1+z),pba->m_ncdm_in_eV[n_ncdm],&qmax),
      pba->error_message,
      pba->error_message);

      // eps = sqrt(q*q*qmax*qmax*(1+z)*(1+z) + pba->m_ncdm_in_eV[n_ncdm]*pba->m_ncdm_in_eV[n_ncdm]);
      // eps = sqrt(q*q*qmax*qmax*(1+z)*(1+z) + pba->m_ncdm_in_eV[n_ncdm]*pba->m_ncdm_in_eV[n_ncdm]);
      // eps = sqrt(q*q*qmax*qmax*(1+z)*(1+z));//we define q in units of qmax.
      // if(mu_ncdm > eps)
      eps = sqrt(q*q*qmax*qmax*(1+z)*(1+z) + pba->m_ncdm_in_eV[n_ncdm]*pba->m_ncdm_in_eV[n_ncdm]); //VP: need to look into this
      // printf("n_ncdm %d eps %e\n", n_ncdm,eps);
      if(pba->ncdm_background_distribution[n_ncdm]==_majoron_){
        *f0=1.0/pow(2*_PI_,3)*(1./(exp((eps-mu_ncdm)/T_ncdm)-1));//bose-einstein

      }
      else{
        *f0=1.0/pow(2*_PI_,3)*(1./(exp((q*qmax*(1+z)-mu_ncdm)/T_ncdm)+1)); //frozen fermi-dirac distribution
      }
      // if(n_ncdm == 1)printf("z %e *f0 %e\n",z,*f0);
      if(*f0 == 0 || isnan(*f0) || *f0 < 0 || *f0 < 1e-40) *f0 = 1e-40; //to avoid bug; eps/T_ncdm can become too big for the exponential when m>>T. we could probably improve that but it works.
      // if(*f0 < 0)printf("n_ncdm %d *f0 %e z %e eps %e mu_ncdm %e T_ncdm %e exp((eps-mu_ncdm)/T_ncdm) %e\n",n_ncdm,*f0,z,eps,mu_ncdm,T_ncdm,exp((eps-mu_ncdm)/T_ncdm));
      // if(1+z<1.5 && n_ncdm == 0)printf("here (1+z) %e ncdm %d mu_ncdm %e Tnu %e q %e eps %e Mncm %e  exp((eps-mu_ncdm)/T_ncdm) %e f0 %e\n",1+z,n_ncdm,mu_ncdm,T_ncdm, q,eps,pba->m_ncdm_in_eV[n_ncdm],exp((eps-mu_ncdm)/T_ncdm),*f0);

      }
    }

  return _SUCCESS_;
}

int background_ncdm_distribution_at_eps(
                                 struct background * pba,
                                 double z,
                                 int n_ncdm,
                                 double eps,
                                 double * f0
                                 ) {

   double T_ncdm,mu_ncdm,p,E;
   //VP: here we define the background distribution
   interpolate_T_and_mu_at_z(pba,n_ncdm,z,&T_ncdm,&mu_ncdm);
   // printf("eps*(1+z)\n", eps*(1+z));
   E = eps*(1+z);
   if(pba->ncdm_background_distribution[n_ncdm]==_majoron_){
     *f0=1.0/pow(2*_PI_,3)*(1./(exp((E-mu_ncdm)/T_ncdm)-1));//bose-einstein

   }
   else{
     p = sqrt(E*E-pba->m_ncdm_in_eV[n_ncdm]*pba->m_ncdm_in_eV[n_ncdm]);
     *f0=1.0/pow(2*_PI_,3)*(1./(exp((p-mu_ncdm)/T_ncdm)+1)); //frozen fermi-dirac distribution
   }
   // if(n_ncdm == 1)printf("z %e *f0 %e\n",z,*f0);
   if(*f0 == 0 || isnan(*f0) || *f0 < 0 || *f0 < 1e-40) *f0 = 1e-40; //to avoid bug; eps/T_ncdm can become too big for the exponential when m>>T. we could probably improve that but it works.
   // if(*f0 == 1e-40)printf("n_ncdm %d *f0 %e z %e eps %e mu_ncdm %e T_ncdm %e exp((eps-mu_ncdm)/T_ncdm) %e\n",n_ncdm,*f0,z,E,mu_ncdm,T_ncdm,exp((E-mu_ncdm)/T_ncdm));
   // if(1+z<1.5 && n_ncdm == 0)printf("here (1+z) %e ncdm %d mu_ncdm %e Tnu %e q %e eps %e Mncm %e  exp((eps-mu_ncdm)/T_ncdm) %e f0 %e\n",1+z,n_ncdm,mu_ncdm,T_ncdm, q,eps,pba->m_ncdm_in_eV[n_ncdm],exp((eps-mu_ncdm)/T_ncdm),*f0);
   return _SUCCESS_;

}
int interpolate_T_and_mu_at_z(struct background *pba,int n_ncdm, double z,double *T_ncdm, double *mu_ncdm){
  //VP:This function is used to interpolate and extrapolate T and mu
  //First: check if spline interpolation is possible
  //If not: extrapolate using analytical formulae
 int last_index;
 int start_indx=20;

  if (z >= pba->z_maj[pba->len_maj] && z <= pba->z_maj[start_indx]){

    if(n_ncdm == pba->entry_is_M_phi){
      /** - interpolate from pre-computed table with array_interpolate() */
      class_call(array_interpolate_spline(
                                          pba->z_maj,
                                          pba->len_maj+1,
                                          pba->T_maj,
                                          pba->ddT_maj,
                                          1,
                                          z,
                                          &last_index,
                                          T_ncdm,
                                          1,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);
      class_call(array_interpolate_spline(
                                          pba->z_maj,
                                          pba->len_maj+1,
                                          pba->Mu_maj,
                                          pba->ddMu_maj,
                                          1,
                                          z,
                                          &last_index,
                                          mu_ncdm,
                                          1,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);

    }else{
      /** - interpolate from pre-computed table with array_interpolate() */
      class_call(array_interpolate_spline(
                                          pba->z_maj,
                                          pba->len_maj+1,
                                          pba->T_nu,
                                          pba->ddT_nu,
                                          1,
                                          z,
                                          &last_index,
                                          T_ncdm,
                                          1,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);
      class_call(array_interpolate_spline(
                                          pba->z_maj,
                                          pba->len_maj+1,
                                          pba->Mu_nu,
                                          pba->ddMu_nu,
                                          1,
                                          z,
                                          &last_index,
                                          mu_ncdm,
                                          1,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);
    }
  }else if(z > pba->z_maj[start_indx]){
        //at early times
        if(n_ncdm == pba->entry_is_M_phi){
          *mu_ncdm = (pba->Mu_maj[start_indx])*(1+z)/(1+pba->z_maj[start_indx]);
          *T_ncdm = pba->T_maj[start_indx]*(1+z)/(1+pba->z_maj[start_indx]);
        }
        else{
          // *mu_ncdm = (pba->Mu_nu[0]-pba->m_ncdm_in_eV[n_ncdm])*(1+z)/(1+pba->z_maj[0]);
          *mu_ncdm = (pba->Mu_nu[start_indx])*(1+z)/(1+pba->z_maj[start_indx]);
          *T_ncdm = pba->T_nu[start_indx]*(1+z)/(1+pba->z_maj[start_indx]);
        }
  }
  else{//at late times
      if(n_ncdm == pba->entry_is_M_phi){
        //assume relativistic:
        *T_ncdm = pba->T_maj[pba->len_maj]*(1+z)/(1+pba->z_maj[pba->len_maj]);
        // *mu_ncdm = (pba->Mu_maj[pba->len_maj]-pba->m_ncdm_in_eV[n_ncdm])*(1+z)/(1+pba->z_maj[pba->len_maj]);
        *mu_ncdm = (pba->Mu_maj[pba->len_maj])*(1+z)/(1+pba->z_maj[pba->len_maj]);
        if(*T_ncdm < 3./20*pba->m_ncdm_in_eV[n_ncdm]){
          if(pba->z_nrel[n_ncdm] <= 0)pba->z_nrel[n_ncdm] = z;
          //check that we were correct.
          // *mu_ncdm = (pba->Mu_maj[pba->len_maj])*(1+z_nrel[n_ncdm])/(1+pba->z_maj[pba->len_maj])*pow((1+z)/(1+z_nrel[n_ncdm]),2);//if we include the mass leads to big discontinuity.//we decided that at late times we don't care about mu because rho is so small.
          *T_ncdm = pba->T_maj[pba->len_maj]*(1+pba->z_nrel[n_ncdm])/(1+pba->z_maj[pba->len_maj])*pow((1+z)/(1+pba->z_nrel[n_ncdm]),2);
        }

      // printf("z %e \n",z,*T_ncdm/pba->m_ncdm_in_eV[n_ncdm]);


    }else{
      //assume relativistic:
      *T_ncdm = pba->T_nu[pba->len_maj]*(1+z)/(1+pba->z_maj[pba->len_maj]);
      *mu_ncdm = (pba->Mu_nu[pba->len_maj])*(1+z)/(1+pba->z_maj[pba->len_maj]);
        if(*T_ncdm < 3./20*pba->m_ncdm_in_eV[n_ncdm]){
          if(pba->z_nrel[n_ncdm] <= 0)pba->z_nrel[n_ncdm] = z;
          //check that we were correct.
          // *mu_ncdm = (pba->Mu_nu[pba->len_maj])*(1+pba->z_nrel[n_ncdm])/(1+pba->z_maj[pba->len_maj])*pow((1+z)/(1+pba->z_nrel[n_ncdm]),2);//if we include the mass leads to big discontinuity; is that ok?
          // *mu_ncdm = (pba->Mu_nu[pba->len_maj]-pba->m_ncdm_in_eV[n_ncdm])*(1+pba->z_nrel[n_ncdm])/(1+pba->z_maj[pba->len_maj])*pow((1+z)/(1+pba->z_nrel[n_ncdm]),2);
          // *T_ncdm = pba->T_nu[pba->len_maj]*(1+pba->z_nrel[n_ncdm])/(1+pba->z_maj[pba->len_maj])*pow((1+z)/(1+pba->z_nrel[n_ncdm]),2);
          //VP: HERE NEED TO UNDERSTAND WHY THE CODE BREAKS IF T AND MU EXTRAPOLATION CHANGES.
        }
    }

    // printf("pba->z_nrel[n_ncdm = %d] %e\n",n_ncdm,pba->z_nrel[n_ncdm]);
  }

  return _SUCCESS_;
}
/**
 * This function is only used for the purpose of finding optimal
 * quadrature weights. The logic is: if we can accurately convolve
 * f0(q) with this function, then we can convolve it accurately with
 * any other relevant function.
 *
 * @param pbadist Input:  structure containing all background parameters
 * @param q       Input:  momentum
 * @param test    Output: value of the test function test(q)
 */

int background_ncdm_test_function(
                                  void * pbadist,
                                  double q,
                                  double * test
                                  ) {

  double c = 2.0/(3.0*_zeta3_);
  double d = 120.0/(7.0*pow(_PI_,4));
  double e = 2.0/(45.0*_zeta5_);

  /** Using a + bq creates problems for otherwise acceptable distributions
      which diverges as \f$ 1/r \f$ or \f$ 1/r^2 \f$ for \f$ r\to 0 \f$*/
  *test = pow(2.0*_PI_,3)/6.0*(c*q*q-d*q*q*q-e*q*q*q*q);

  return _SUCCESS_;
}

/**
 * This function finds optimal quadrature weights for each ncdm
 * species
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */


int get_q_max(struct background *pba, int n_ncdm, double a, double M,double * qmax){
  //VP: extract the maximum comoving momentum at a
  //minimum is always 0.
  double T_ncdm;
  double mu_ncdm;
  class_call(interpolate_T_and_mu_at_z(pba,n_ncdm,1./a-1.,&T_ncdm,&mu_ncdm),
  pba->error_message,
  pba->error_message);
  if(20*T_ncdm > 3 *M){
    *qmax = pow(20*20*T_ncdm*T_ncdm-(M*M),0.5)*a;
  }else{
    // *qmax = pow(3*3*M*M-(M*M),0.5)*a;
    *qmax = 20*T_ncdm*a; //VP: need to look into this
  }
 *qmax = 20*T_ncdm*a;
  return _SUCCESS_;
}


int interpolate_background_ncdm_distribution(struct background *pba, int n_ncdm, double *qtable,double qsize, double z, double *ftable) {
  //VP: this function fills a table of f_ncdm(q) at a.
  //the qtable is given in unit of qmax.
  //qmax is then computed within background_ncdm_distribution.

  int index_q;
  struct background_parameters_for_distributions pbadist;
  pbadist.pba = pba;
  pbadist.n_ncdm = n_ncdm;
  double f0;
  for(index_q = 0; index_q < qsize; index_q++){

    class_call(background_ncdm_distribution(
                                     &pbadist,
                                     qtable[index_q],
                                     &f0,
                                     z),
                 pba->error_message,
                 pba->error_message);
    ftable[index_q]= f0;
  }

  return _SUCCESS_;

}

// SJW -- ADD HERE
int background_ncdm_init(
                         struct precision *ppr,
                         struct background *pba
                         ) {

  int index_q, k,tolexp,row,status,filenum;
  double f0m2,f0m1,f0,f0p1,f0p2,dq,q,df0dq,tmp1,tmp2;
  struct background_parameters_for_distributions pbadist;
  FILE *psdfile;

  pbadist.pba = pba;

  /* Allocate pointer arrays: */
  class_alloc(pba->q_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->w_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->dlnf0_dlnq_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->f_ncdm_bg, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->f_ncdm, sizeof(double*)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->z_nrel, sizeof(double*)*pba->N_ncdm,pba->error_message);

  /* Allocate pointers: */
  class_alloc(pba->q_size_ncdm,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->q_size_ncdm_bg,sizeof(int)*pba->N_ncdm,pba->error_message);
  class_alloc(pba->factor_ncdm,sizeof(double)*pba->N_ncdm,pba->error_message);



   /// SJW -- Implimenting background functions....
   int test, lenIndx;
   //printf("At function... \n");
   if(pba->M_phi > 0){
     //we do this if we have a majoron
     class_alloc(pba->z_maj, sizeof(double*)*1000,pba->error_message);
     class_alloc(pba->T_maj, sizeof(double*)*1000,pba->error_message);
     class_alloc(pba->T_nu, sizeof(double*)*1000,pba->error_message);
     class_alloc(pba->Mu_maj, sizeof(double*)*1000,pba->error_message);
     class_alloc(pba->Mu_nu, sizeof(double*)*1000,pba->error_message);
     //VP: add these tables for purpose of spline interpolation
     class_alloc(pba->ddT_maj, sizeof(double*)*1000,pba->error_message);
     class_alloc(pba->ddT_nu, sizeof(double*)*1000,pba->error_message);
     class_alloc(pba->ddMu_maj, sizeof(double*)*1000,pba->error_message);
     class_alloc(pba->ddMu_nu, sizeof(double*)*1000,pba->error_message);


     test = background_MB_approx(pba, &lenIndx);
     pba->len_maj = lenIndx;

     //VP: once we know the size of the table, we can "realloc"; i.e. removes useless space in tables.
     pba->z_maj = realloc(pba->z_maj, sizeof(double*)*(pba->len_maj+1));
     pba->T_maj = realloc(pba->T_maj, sizeof(double*)*(pba->len_maj+1));
     pba->T_nu = realloc(pba->T_nu, sizeof(double*)*(pba->len_maj+1));
     pba->Mu_maj = realloc(pba->Mu_maj, sizeof(double*)*(pba->len_maj+1));
     pba->Mu_nu = realloc(pba->Mu_nu, sizeof(double*)*(pba->len_maj+1));
     pba->ddT_maj = realloc(pba->ddT_maj, sizeof(double*)*(pba->len_maj+1));
     pba->ddT_nu = realloc(pba->ddT_nu, sizeof(double*)*(pba->len_maj+1));
     pba->ddMu_maj = realloc(pba->ddMu_maj, sizeof(double*)*(pba->len_maj+1));
     pba->ddMu_nu = realloc(pba->ddMu_nu, sizeof(double*)*(pba->len_maj+1));


     //VP: getting ready to interpolate
     class_call(array_spline_table_lines(pba->z_maj,
                                         pba->len_maj+1,
                                         pba->T_maj,
                                         1,
                                         pba->ddT_maj,
                                         _SPLINE_EST_DERIV_,
                                         pba->error_message),
      pba->error_message,
      pba->error_message);
     class_call(array_spline_table_lines(pba->z_maj,
                                         pba->len_maj+1,
                                         pba->Mu_maj,
                                         1,
                                         pba->ddMu_maj,
                                         _SPLINE_EST_DERIV_,
                                         pba->error_message),
      pba->error_message,
      pba->error_message);
     class_call(array_spline_table_lines(pba->z_maj,
                                         pba->len_maj+1,
                                         pba->Mu_nu,
                                         1,
                                         pba->ddMu_nu,
                                         _SPLINE_EST_DERIV_,
                                         pba->error_message),
      pba->error_message,
      pba->error_message);
     class_call(array_spline_table_lines(pba->z_maj,
                                         pba->len_maj+1,
                                         pba->T_nu,
                                         1,
                                         pba->ddT_nu,
                                         _SPLINE_EST_DERIV_,
                                         pba->error_message),
      pba->error_message,
      pba->error_message);

     //
     //
     // for (int ii=0; ii <= pba->len_maj; ii++){
     //     printf("HERE: %d \t %e \t %e \t %e \t %e \t %e \n", ii, pba->z_maj[ii], pba->T_maj[ii], pba->T_nu[ii], pba->Mu_maj[ii], pba->Mu_nu[ii]);
     // }
     // exit(0);
  //
   }


  for(k=0, filenum=0; k<pba->N_ncdm; k++){
    pbadist.n_ncdm = k;
    pbadist.q = NULL;
    pbadist.tablesize = 0;
    /*Do we need to read in a file to interpolate the distribution function? */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      psdfile = fopen(pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_,"r");
      class_test(psdfile == NULL,pba->error_message,
                 "Could not open file %s!",pba->ncdm_psd_files+filenum*_ARGUMENT_LENGTH_MAX_);
      // Find size of table:
      for (row=0,status=2; status==2; row++){
        status = fscanf(psdfile,"%lf %lf",&tmp1,&tmp2);
      }
      rewind(psdfile);
      pbadist.tablesize = row-1;

      /*Allocate room for interpolation table: */
      class_alloc(pbadist.q,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      class_alloc(pbadist.d2f0,sizeof(double)*pbadist.tablesize,pba->error_message);
      for (row=0; row<pbadist.tablesize; row++){
        status = fscanf(psdfile,"%lf %lf",
                        &pbadist.q[row],&pbadist.f0[row]);
        //		printf("(q,f0) = (%g,%g)\n",pbadist.q[row],pbadist.f0[row]);
      }
      fclose(psdfile);
      /* Call spline interpolation: */
      class_call(array_spline_table_lines(pbadist.q,
                                          pbadist.tablesize,
                                          pbadist.f0,
                                          1,
                                          pbadist.d2f0,
                                          _SPLINE_EST_DERIV_,
                                          pba->error_message),
                 pba->error_message,
                 pba->error_message);
      filenum++;
    }

    /* Handle perturbation qsampling: */
    if (pba->ncdm_quadrature_strategy[k]==qm_auto){
      /** Automatic q-sampling for this species */
      class_alloc(pba->q_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm[k],_QUADRATURE_MAX_*sizeof(double),pba->error_message);

      class_call(get_qsampling(pba->q_ncdm[k],
                               pba->w_ncdm[k],
                               &(pba->q_size_ncdm[k]),
                               _QUADRATURE_MAX_,
                               ppr->tol_ncdm,
                               pbadist.q,
                               pbadist.tablesize,
                               background_ncdm_test_function,
                               background_ncdm_distribution,
                               &pbadist,
                               pba->error_message),
                 pba->error_message,
                 pba->error_message);
      pba->q_ncdm[k]=realloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));
      pba->w_ncdm[k]=realloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double));


      if (pba->background_verbose > 0)
        printf("ncdm species i=%d sampled with %d points for purpose of perturbation integration\n",
               k+1,
               pba->q_size_ncdm[k]);

      /* Handle background q_sampling: */
      class_alloc(pba->q_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm_bg[k],_QUADRATURE_MAX_BG_*sizeof(double),pba->error_message);

      // printf("ppr->tol_ncdm_bg %e\n", ppr->tol_ncdm_bg);
      class_call(get_qsampling(pba->q_ncdm_bg[k],
                               pba->w_ncdm_bg[k],
                               &(pba->q_size_ncdm_bg[k]),
                               _QUADRATURE_MAX_BG_,
                               ppr->tol_ncdm_bg,
                               pbadist.q,
                               pbadist.tablesize,
                               background_ncdm_test_function,
                               background_ncdm_distribution,
                               &pbadist,
                               pba->error_message),
                 pba->error_message,
                 pba->error_message);


      pba->q_ncdm_bg[k]=realloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));
      pba->w_ncdm_bg[k]=realloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double));

      /** - in verbose mode, inform user of number of sampled momenta
          for background quantities */
      if (pba->background_verbose > 0)
        printf("ncdm species i=%d sampled with %d points for purpose of background integration\n",
               k+1,
               pba->q_size_ncdm_bg[k]);
    }
    else{
      //VP: for majoron and our neutrinos, we will always be here.
      /** Manual q-sampling for this species. Same sampling used for both perturbation and background sampling, since this will usually be a high precision setting anyway */
      pba->q_size_ncdm_bg[k] = pba->ncdm_input_q_size_bg[k];
      pba->q_size_ncdm[k] = pba->ncdm_input_q_size[k];
      class_alloc(pba->q_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm_bg[k],pba->q_size_ncdm_bg[k]*sizeof(double),pba->error_message);
      class_alloc(pba->q_ncdm[k],pba->q_size_ncdm[k]*sizeof(double),pba->error_message);
      class_alloc(pba->w_ncdm[k],pba->q_size_ncdm[k]*sizeof(double),pba->error_message);
      class_call(get_qsampling_manual(pba->q_ncdm[k],
                                      pba->w_ncdm[k],
                                      pba->q_size_ncdm[k],
                                      pba->ncdm_qmax[k],
                                      pba->ncdm_quadrature_strategy[k],
                                      pbadist.q,
                                      pbadist.tablesize,
                                      background_ncdm_distribution,
                                      &pbadist,
                                      pba->error_message),
                 pba->error_message,
                 pba->error_message);
      class_call(get_qsampling_manual(pba->q_ncdm_bg[k],
                                      pba->w_ncdm_bg[k],
                                      pba->q_size_ncdm_bg[k],
                                      pba->ncdm_qmax[k],
                                      pba->ncdm_quadrature_strategy[k],
                                      pbadist.q,
                                      pbadist.tablesize,
                                      background_ncdm_distribution,
                                      &pbadist,
                                      pba->error_message),
                 pba->error_message,
                 pba->error_message);

      /** - in verbose mode, inform user of number of sampled momenta
          for background quantities */
      if (pba->background_verbose > 0)
        printf("ncdm species i=%d sampled with %d points for purpose of background andperturbation integration using the manual method\n",
               k+1,
               pba->q_size_ncdm[k]);
    }

    class_alloc(pba->dlnf0_dlnq_ncdm[k],
                pba->q_size_ncdm[k]*sizeof(double),
                pba->error_message);
    class_alloc(pba->f_ncdm_bg[k],
                pba->q_size_ncdm_bg[k]*sizeof(double),
                pba->error_message);
    class_alloc(pba->f_ncdm[k],
                pba->q_size_ncdm[k]*sizeof(double),
                pba->error_message);


    if(pba->ncdm_background_distribution[k] == _fermi_dirac_){
    for (index_q=0; index_q<pba->q_size_ncdm[k]; index_q++) {
      q = pba->q_ncdm[k][index_q];
      class_call(background_ncdm_distribution(&pbadist,q,&f0,0),
                 pba->error_message,pba->error_message);

      //we correct the integration weights if we deal with majorons
      //Loop to find appropriate dq:
      for(tolexp=_PSD_DERIVATIVE_EXP_MIN_; tolexp<_PSD_DERIVATIVE_EXP_MAX_; tolexp++){

        if (index_q == 0){
          dq = MIN((0.5-ppr->smallest_allowed_variation)*q,2*exp(tolexp)*(pba->q_ncdm[k][index_q+1]-q));
        }
        else if (index_q == pba->q_size_ncdm[k]-1){
          dq = exp(tolexp)*2.0*(pba->q_ncdm[k][index_q]-pba->q_ncdm[k][index_q-1]);
        }
        else{
          dq = exp(tolexp)*(pba->q_ncdm[k][index_q+1]-pba->q_ncdm[k][index_q-1]);
        }

        class_call(background_ncdm_distribution(&pbadist,q-2*dq,&f0m2,0),
                   pba->error_message,pba->error_message);
        class_call(background_ncdm_distribution(&pbadist,q+2*dq,&f0p2,0),
                   pba->error_message,pba->error_message);

        if (fabs((f0p2-f0m2)/f0)>sqrt(ppr->smallest_allowed_variation)) break;
      }

      class_call(background_ncdm_distribution(&pbadist,q-dq,&f0m1,0),
                 pba->error_message,pba->error_message);
      class_call(background_ncdm_distribution(&pbadist,q+dq,&f0p1,0),
                 pba->error_message,pba->error_message);
      //5 point estimate of the derivative:
      df0dq = (+f0m2-8*f0m1+8*f0p1-f0p2)/12.0/dq;
      //printf("df0dq[%g] = %g. dlf=%g ?= %g. f0 =%g.\n",q,df0dq,q/f0*df0dq,
      //Avoid underflow in extreme tail:
      if (fabs(f0)==0.)
        pba->dlnf0_dlnq_ncdm[k][index_q] = -q; /* valid for whatever f0 with exponential tail in exp(-q) */
      else
        pba->dlnf0_dlnq_ncdm[k][index_q] = q/f0*df0dq;
      // printf("index_q %d q/f0*df0dq %e\n", index_q,q/f0*df0dq);
      }





      pba->factor_ncdm[k]=pba->deg_ncdm[k]*4*_PI_*pow(pba->T_cmb*pba->T_ncdm[k]*_k_B_,4)*8*_PI_*_G_
        /3./pow(_h_P_/2./_PI_,3)/pow(_c_,7)*_Mpc_over_m_*_Mpc_over_m_;//energy is in unit of T_ncdm



    }else if(pba->ncdm_background_distribution[k] == _fermi_dirac_v2_ || pba->ncdm_background_distribution[k] == _majoron_){
      //VP: we will compute dlnf0 later, as a function of time.
      pba->factor_ncdm[k]=pba->deg_ncdm[k]*4*_PI_*pow(_eV_,4)*8*_PI_*_G_
        /3./pow(_h_P_/2./_PI_,3)/pow(_c_,7)*_Mpc_over_m_*_Mpc_over_m_; //in this case, units are eV. we convert to Mpc.
    }

    pba->z_nrel[k] = -1;//initialization; will be attributed later.
    /* If allocated, deallocate interpolation table:  */
    if ((pba->got_files!=NULL)&&(pba->got_files[k]==_TRUE_)){
      free(pbadist.q);
      free(pbadist.f0);
      free(pbadist.d2f0);
    }
  }


  return _SUCCESS_;
}

/**
 * For a given ncdm species: given the quadrature weights, the mass
 * and the redshift, find background quantities by a quick weighted
 * sum over.  Input parameters passed as NULL pointers are not
 * evaluated for speed-up
 *
 * @param qvec     Input: sampled momenta
 * @param wvec     Input: quadrature weights
 * @param qsize    Input: number of momenta/weights
 * @param M        Input: mass
 * @param factor   Input: normalization factor for the p.s.d.
 * @param z        Input: redshift
 * @param n        Output: number density
 * @param rho      Output: energy density
 * @param p        Output: pressure
 * @param drho_dM  Output: derivative used in next function
 * @param pseudo_p Output: pseudo-pressure used in perturbation module for fluid approx
 *
 */

int background_ncdm_momenta(
                            /* Only calculate for non-NULL pointers: */
                            double * qvec,
                            double * wvec,
                            double *fvec,
                            int qsize,
                            double M,
                            double qmax,
                            double factor,
                            double z,
                            int n_ncdm,
                            double * n,
                            double * rho, // density
                            double * p,   // pressure
                            double * drho_dM,  // d rho / d M used in next function
                            double * pseudo_p  // pseudo-p used in ncdm fluid approx
                            ) {

  int index_q;
  double epsilon;
  double q2;
  double factor2;
  /** Summary: */

  /** - rescale normalization at given redshift */
  factor2 = factor*pow(1+z,4);

  /** - initialize quantities */
  if (n!=NULL) *n = 0.;
  if (rho!=NULL) *rho = 0.;
  if (p!=NULL) *p = 0.;
  if (drho_dM!=NULL) *drho_dM = 0.;
  if (pseudo_p!=NULL) *pseudo_p = 0.;

  /** - loop over momenta */
  for (index_q=0; index_q<qsize; index_q++) {

    /* squared momentum */
    //VP: q in units of qmax; and so is wvec (which is dq/qmax).
    q2 = qvec[index_q]*qvec[index_q]*qmax*qmax;

    epsilon = sqrt(q2+M*M/(1.+z)/(1.+z));

    /* integrand of the various quantities */
    if (n!=NULL) *n += q2*wvec[index_q]*qmax*fvec[index_q];
    if (rho!=NULL) *rho += q2*epsilon*wvec[index_q]*qmax*fvec[index_q];
    if (p!=NULL) *p += q2*q2/3./epsilon*wvec[index_q]*qmax*fvec[index_q];
    if (drho_dM!=NULL) *drho_dM += q2*M/(1.+z)/(1.+z)/epsilon*wvec[index_q]*qmax*fvec[index_q];
    if (pseudo_p!=NULL) *pseudo_p += pow(q2/epsilon,3)/3.0*wvec[index_q]*qmax*fvec[index_q];
  }

  /** - adjust normalization */
  if (n!=NULL) *n *= factor2/(1.+z);
  if (rho!=NULL) {
    *rho *= factor2;
  }
  if (p!=NULL) {
    *p *= factor2;
  }
  if (drho_dM!=NULL) *drho_dM *= factor2;
  if (pseudo_p!=NULL) *pseudo_p *=factor2;

  return _SUCCESS_;
}

/**
 * When the user passed the density fraction Omega_ncdm or
 * omega_ncdm in input but not the mass, infer the mass with Newton iteration method.
 *
 * @param ppr    Input: precision structure
 * @param pba    Input/Output: background structure
 * @param n_ncdm Input: index of ncdm species
 */

int background_ncdm_M_from_Omega(
                                 struct precision *ppr,
                                 struct background *pba,
                                 int n_ncdm
                                 ) {
  double rho0,rho,n,M,deltaM,drhodM, T_ncdm, mu_ncdm,qmax;
  int iter,maxiter=50;
  int index_q;

  rho0 = pba->H0*pba->H0*pba->Omega0_ncdm[n_ncdm]; /*Remember that rho is defined such that H^2=sum(rho_i) */
  M = 0.0;

  if(pba->ncdm_background_distribution[n_ncdm] == _fermi_dirac_v2_ || pba->ncdm_background_distribution[n_ncdm] == _majoron_){
    class_call(interpolate_background_ncdm_distribution(pba,n_ncdm,pba->q_ncdm_bg[n_ncdm],pba->q_size_ncdm_bg[n_ncdm],0.0,pba->f_ncdm_bg[n_ncdm]),
    pba->error_message,
    pba->error_message);
    class_call(get_q_max(pba,n_ncdm,1.,M,&qmax),
    pba->error_message,
    pba->error_message);
  }else{
    for(index_q = 0; index_q < pba->q_size_ncdm_bg[n_ncdm]; index_q++){
    pba->f_ncdm_bg[n_ncdm][index_q] = 1;//f_ncdm is already included in w_ncdm_bg in the case of standard neutrinos.
    }
    qmax=1;
  }


  background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                          pba->w_ncdm_bg[n_ncdm],
                          pba->f_ncdm_bg[n_ncdm],
                          pba->q_size_ncdm_bg[n_ncdm],
                          M,
                          qmax,
                          pba->factor_ncdm[n_ncdm],
                          0.,
                          n_ncdm,
                          &n,
                          &rho,
                          NULL,
                          NULL,
                          NULL);

  /* Is the value of Omega less than a massless species?*/
  class_test(rho0<rho,pba->error_message,
             "The value of Omega for the %dth species, %g, is less than for a massless species! It should be atleast %g. Check your input.",
             n_ncdm,pba->Omega0_ncdm[n_ncdm],pba->Omega0_ncdm[n_ncdm]*rho/rho0);

  /* In the strict NR limit we have rho = n*(M) today, giving a zeroth order guess: */
  M = rho0/n; /* This is our guess for M. */
  for (iter=1; iter<=maxiter; iter++){

    if(pba->ncdm_background_distribution[n_ncdm] == _fermi_dirac_v2_ || pba->ncdm_background_distribution[n_ncdm] == _majoron_){
      class_call(interpolate_background_ncdm_distribution(pba,n_ncdm,pba->q_ncdm_bg[n_ncdm],pba->q_size_ncdm_bg[n_ncdm],0.0,pba->f_ncdm_bg[n_ncdm]),
      pba->error_message,
      pba->error_message);
      class_call(get_q_max(pba,n_ncdm,1.,M,&qmax),
      pba->error_message,
      pba->error_message);
      //need to check units of M.
    }else{
      for(index_q = 0; index_q < pba->q_size_ncdm_bg[n_ncdm]; index_q++){
      pba->f_ncdm_bg[n_ncdm][index_q] = 1;//f_ncdm is already included in w_ncdm_bg in the case of standard neutrinos.
      }
      qmax =1;
    }


    /* Newton iteration. First get relevant quantities at M: */
    background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                            pba->w_ncdm_bg[n_ncdm],
                            pba->f_ncdm_bg[n_ncdm],
                            pba->q_size_ncdm_bg[n_ncdm],
                            M,
                            qmax,
                            pba->factor_ncdm[n_ncdm],
                            0.,
                            n_ncdm,
                            NULL,
                            &rho,
                            NULL,
                            &drhodM,
                            NULL);

    deltaM = (rho0-rho)/drhodM; /* By definition of the derivative */
    if ((M+deltaM)<0.0) deltaM = -M/2.0; /* Avoid overshooting to negative M value. */
    M += deltaM; /* Update value of M.. */
    if (fabs(deltaM/M)<ppr->tol_M_ncdm){
      /* Accuracy reached.. */
      pba->M_ncdm[n_ncdm] = M;
      break;
    }
  }
  class_test(iter>=maxiter,pba->error_message,
             "Newton iteration could not converge on a mass for some reason.");
  return _SUCCESS_;
}

/**
 *  This function integrates the background over time, allocates and
 *  fills the background table
 *
 * @param ppr Input: precision structure
 * @param pba Input/Output: background structure
 */

int background_solve(
                     struct precision *ppr,
                     struct background *pba
                     ) {

  /** Summary: */

  /** - define local variables */

  /* contains all quantities relevant for the integration algorithm */
  struct generic_integrator_workspace gi;
  /* parameters and workspace for the background_derivs function */
  struct background_parameters_and_workspace bpaw;
  /* a growing table (since the number of time steps is not known a priori) */
  growTable gTable;
  /* needed for growing table */
  double * pData;
  /* needed for growing table */
  void * memcopy_result;
  /* initial conformal time */
  double tau_start;
  /* final conformal time */
  double tau_end;
  /* an index running over bi indices */
  int i;
  /* vector of quantities to be integrated */
  double * pvecback_integration;
  /* vector of all background quantities */
  double * pvecback;
  /* necessary for calling array_interpolate(), but never used */
  int last_index=0, n_ncdm;
  /* comoving radius coordinate in Mpc (equal to conformal distance in flat case) */
  double comoving_radius=0.;

  int lenIndx;

  bpaw.pba = pba;
  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);
  bpaw.pvecback = pvecback;

  /** - allocate vector of quantities to be integrated */
  class_alloc(pvecback_integration,pba->bi_size*sizeof(double),pba->error_message);

  /** - initialize generic integrator with initialize_generic_integrator() */

  /* Size of vector to integrate is (pba->bi_size-1) rather than
   * (pba->bi_size), since tau is not integrated.
   */
  class_call(initialize_generic_integrator((pba->bi_size-1),&gi),
             gi.error_message,
             pba->error_message);

  /** - impose initial conditions with background_initial_conditions() */
  class_call(background_initial_conditions(ppr,pba,pvecback,pvecback_integration),
             pba->error_message,
             pba->error_message);

  /* here tau_end is in fact the initial time (in the next loop
     tau_start = tau_end) */
  tau_end=pvecback_integration[pba->index_bi_tau];

  /** - create a growTable with gt_init() */
  class_call(gt_init(&gTable),
             gTable.error_message,
             pba->error_message);

  /* initialize the counter for the number of steps */
  pba->bt_size=0;

  /** - loop over integration steps: call background_functions(), find step size, save data in growTable with gt_add(), perform one step with generic_integrator(), store new value of tau */

  while (pvecback_integration[pba->index_bi_a] < pba->a_today) {

    tau_start = tau_end;

    /* -> find step size (trying to adjust the last step as close as possible to the one needed to reach a=a_today; need not be exact, difference corrected later) */
    class_call(background_functions(pba,pvecback_integration, pba->short_info, pvecback),
               pba->error_message,
               pba->error_message);

    if ((pvecback_integration[pba->index_bi_a]*(1.+ppr->back_integration_stepsize)) < pba->a_today) {
      tau_end = tau_start + ppr->back_integration_stepsize / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }
    else {
      tau_end = tau_start + (pba->a_today/pvecback_integration[pba->index_bi_a]-1.) / (pvecback_integration[pba->index_bi_a]*pvecback[pba->index_bg_H]);
      /* no possible segmentation fault here: non-zeroness of "a" has been checked in background_functions() */
    }

    class_test((tau_end-tau_start)/tau_start < ppr->smallest_allowed_variation,
               pba->error_message,
               "integration step: relative change in time =%e < machine precision : leads either to numerical error or infinite loop",(tau_end-tau_start)/tau_start);

    /* -> save data in growTable */
    class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
               gTable.error_message,
               pba->error_message);
    pba->bt_size++;

    /* -> perform one step */
    class_call(generic_integrator(background_derivs,
                                  tau_start,
                                  tau_end,
                                  pvecback_integration,
                                  &bpaw,
                                  ppr->tol_background_integration,
                                  ppr->smallest_allowed_variation,
                                  &gi),
               gi.error_message,
               pba->error_message);

    /* -> store value of tau */
    pvecback_integration[pba->index_bi_tau]=tau_end;

  }

  /** - save last data in growTable with gt_add() */
  class_call(gt_add(&gTable,_GT_END_,(void *) pvecback_integration,sizeof(double)*pba->bi_size),
             gTable.error_message,
             pba->error_message);
  pba->bt_size++;


  /* integration finished */

  /** - clean up generic integrator with cleanup_generic_integrator() */
  class_call(cleanup_generic_integrator(&gi),
             gi.error_message,
             pba->error_message);

  /** - retrieve data stored in the growTable with gt_getPtr() */
  class_call(gt_getPtr(&gTable,(void**)&pData),
             gTable.error_message,
             pba->error_message);

  /** - interpolate to get quantities precisely today with array_interpolate() */
  class_call(array_interpolate(
                               pData,
                               pba->bi_size,
                               pba->bt_size,
                               pba->index_bi_a,
                               pba->a_today,
                               &last_index,
                               pvecback_integration,
                               pba->bi_size,
                               pba->error_message),
             pba->error_message,
             pba->error_message);

  /* substitute last line with quantities today */
  for (i=0; i<pba->bi_size; i++)
    pData[(pba->bt_size-1)*pba->bi_size+i]=pvecback_integration[i];

  /** - deduce age of the Universe */
  /* -> age in Gyears */
  pba->age = pvecback_integration[pba->index_bi_time]/_Gyr_over_Mpc_;
  /* -> conformal age in Mpc */
  pba->conformal_age = pvecback_integration[pba->index_bi_tau];
  /* -> contribution of decaying dark matter and dark radiation to the critical density today: */
  if (pba->has_dcdm == _TRUE_){
    pba->Omega0_dcdm = pvecback_integration[pba->index_bi_rho_dcdm]/pba->H0/pba->H0;
  }
  if (pba->has_dr == _TRUE_){
    pba->Omega0_dr = pvecback_integration[pba->index_bi_rho_dr]/pba->H0/pba->H0;
  }

  /** - allocate background tables */
  class_alloc(pba->tau_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->z_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2tau_dz2_table,pba->bt_size * sizeof(double),pba->error_message);

  class_alloc(pba->background_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  class_alloc(pba->d2background_dtau2_table,pba->bt_size * pba->bg_size * sizeof(double),pba->error_message);

  /** - In a loop over lines, fill background table using the result of the integration plus background_functions() */

  for (i=0; i < pba->bt_size; i++) {

    /* -> establish correspondence between the integrated variable and the bg variables */

    pba->tau_table[i] = pData[i*pba->bi_size+pba->index_bi_tau];

    class_test(pData[i*pba->bi_size+pba->index_bi_a] <= 0.,
               pba->error_message,
               "a = %e instead of strictly positiv",pData[i*pba->bi_size+pba->index_bi_a]);

    pba->z_table[i] = pba->a_today/pData[i*pba->bi_size+pba->index_bi_a]-1.;

    pvecback[pba->index_bg_time] = pData[i*pba->bi_size+pba->index_bi_time];
    pvecback[pba->index_bg_conf_distance] = pba->conformal_age - pData[i*pba->bi_size+pba->index_bi_tau];

    if (pba->sgnK == 0) comoving_radius = pvecback[pba->index_bg_conf_distance];
    else if (pba->sgnK == 1) comoving_radius = sin(sqrt(pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(pba->K);
    else if (pba->sgnK == -1) comoving_radius = sinh(sqrt(-pba->K)*pvecback[pba->index_bg_conf_distance])/sqrt(-pba->K);

    pvecback[pba->index_bg_ang_distance] = pba->a_today*comoving_radius/(1.+pba->z_table[i]);
    pvecback[pba->index_bg_lum_distance] = pba->a_today*comoving_radius*(1.+pba->z_table[i]);
    pvecback[pba->index_bg_rs] = pData[i*pba->bi_size+pba->index_bi_rs];

    /* -> compute all other quantities depending only on {B} variables.
       The value of {B} variables in pData are also copied to pvecback.*/



    class_call(background_functions(pba,pData+i*pba->bi_size, pba->long_info, pvecback),
               pba->error_message,
               pba->error_message);

    /* -> compute growth functions (valid in dust universe) */

    /* Normalise D(z=0)=1 and construct f = D_prime/(aHD) */
    pvecback[pba->index_bg_D] = pData[i*pba->bi_size+pba->index_bi_D]/pData[(pba->bt_size-1)*pba->bi_size+pba->index_bi_D];
    pvecback[pba->index_bg_f] = pData[i*pba->bi_size+pba->index_bi_D_prime]/
      (pData[i*pba->bi_size+pba->index_bi_D]*pvecback[pba->index_bg_a]*pvecback[pba->index_bg_H]);

    /* -> write in the table */
    memcopy_result = memcpy(pba->background_table + i*pba->bg_size,pvecback,pba->bg_size*sizeof(double));

    class_test(memcopy_result != pba->background_table + i*pba->bg_size,
               pba->error_message,
               "cannot copy data back to pba->background_table");

  }


  /** - free the growTable with gt_free() */

  class_call(gt_free(&gTable),
             gTable.error_message,
             pba->error_message);

   if(pba->has_ncdm){
     for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {

       // if(pba->ncdm_background_distribution[n_ncdm] == _fermi_dirac_v2_ || pba->ncdm_background_distribution[n_ncdm] == _majoron_){


       /** - ---> second derivative with respect to tau of cb2 */
       class_call(array_spline_table_line_to_line(pba->tau_table,
                                                  pba->bt_size,
                                                  pba->background_table,
                                                  pba->bg_size,
                                                  pba->index_bg_T_ncdm1+n_ncdm,
                                                  pba->index_bg_ddT_ncdm1+n_ncdm,
                                                  _SPLINE_EST_DERIV_,
                                                  pba->error_message),
                  pba->error_message,
                  pba->error_message);


       /** - ---> first derivative with respect to tau of cb2 (using spline interpolation) */
       class_call(array_derive_spline_table_line_to_line(pba->tau_table,
                                                         pba->bt_size,
                                                         pba->background_table,
                                                         pba->bg_size,
                                                         pba->index_bg_T_ncdm1+n_ncdm,
                                                         pba->index_bg_ddT_ncdm1+n_ncdm,
                                                         pba->index_bg_dT_ncdm1+n_ncdm,
                                                         pba->error_message),
                  pba->error_message,
                  pba->error_message);
       /** - ---> second derivative with respect to tau of cb2 */
       class_call(array_spline_table_line_to_line(pba->tau_table,
                                                  pba->bt_size,
                                                  pba->background_table,
                                                  pba->bg_size,
                                                  pba->index_bg_Mu_ncdm1+n_ncdm,
                                                  pba->index_bg_ddMu_ncdm1+n_ncdm,
                                                  _SPLINE_EST_DERIV_,
                                                  pba->error_message),
                  pba->error_message,
                  pba->error_message);


       /** - ---> first derivative with respect to tau of cb2 (using spline interpolation) */
       class_call(array_derive_spline_table_line_to_line(pba->tau_table,
                                                         pba->bt_size,
                                                         pba->background_table,
                                                         pba->bg_size,
                                                         pba->index_bg_Mu_ncdm1+n_ncdm,
                                                         pba->index_bg_ddMu_ncdm1+n_ncdm,
                                                         pba->index_bg_dMu_ncdm1+n_ncdm,
                                                         pba->error_message),
                  pba->error_message,
                  pba->error_message);
     // }
    }
   }



  /** - fill tables of second derivatives (in view of spline interpolation) */
  class_call(array_spline_table_lines(pba->z_table,
                                      pba->bt_size,
                                      pba->tau_table,
                                      1,
                                      pba->d2tau_dz2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);

  class_call(array_spline_table_lines(pba->tau_table,
                                      pba->bt_size,
                                      pba->background_table,
                                      pba->bg_size,
                                      pba->d2background_dtau2_table,
                                      _SPLINE_EST_DERIV_,
                                      pba->error_message),
             pba->error_message,
             pba->error_message);






  /** - compute remaining "related parameters" */

  /**  - so-called "effective neutrino number", computed at earliest
      time in interpolation table. This should be seen as a
      definition: Neff is the equivalent number of
      instantaneously-decoupled neutrinos accounting for the
      radiation density, beyond photons */

  pba->Neff = (pba->background_table[pba->index_bg_Omega_r]
               *pba->background_table[pba->index_bg_rho_crit]
               -pba->background_table[pba->index_bg_rho_g])
    /(7./8.*pow(4./11.,4./3.)*pba->background_table[pba->index_bg_rho_g]);

  /** - done */
  if (pba->background_verbose > 0) {
    printf(" -> age = %f Gyr\n",pba->age);
    printf(" -> conformal age = %f Mpc\n",pba->conformal_age);
  }

  if (pba->background_verbose > 2) {
    printf(" -> pba->Neff = %f\n",pba->Neff);
    if ((pba->has_dcdm == _TRUE_)&&(pba->has_dr == _TRUE_)){
      printf("    Decaying Cold Dark Matter details: (DCDM --> DR)\n");
      printf("     -> Omega0_dcdm = %f\n",pba->Omega0_dcdm);
      printf("     -> Omega0_dr = %f\n",pba->Omega0_dr);
      printf("     -> Omega0_dr+Omega0_dcdm = %f, input value = %f\n",
             pba->Omega0_dr+pba->Omega0_dcdm,pba->Omega0_dcdmdr);
      printf("     -> Omega_ini_dcdm/Omega_b = %f\n",pba->Omega_ini_dcdm/pba->Omega0_b);
    }
    if (pba->has_scf == _TRUE_){
      printf("    Scalar field details:\n");
      printf("     -> Omega_scf = %g, wished %g\n",
             pvecback[pba->index_bg_rho_scf]/pvecback[pba->index_bg_rho_crit], pba->Omega0_scf);
      if(pba->has_lambda == _TRUE_)
        printf("     -> Omega_Lambda = %g, wished %g\n",
               pvecback[pba->index_bg_rho_lambda]/pvecback[pba->index_bg_rho_crit], pba->Omega0_lambda);
      printf("     -> parameters: [lambda, alpha, A, B] = \n");
      printf("                    [");
      for (i=0; i<pba->scf_parameters_size-1; i++){
        printf("%.3f, ",pba->scf_parameters[i]);
      }
      printf("%.3f]\n",pba->scf_parameters[pba->scf_parameters_size-1]);
    }
  }

  /**  - total matter, radiation, dark energy today */
  pba->Omega0_m = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_Omega_m];
  pba->Omega0_r = pba->background_table[(pba->bt_size-1)*pba->bg_size+pba->index_bg_Omega_r];
  pba->Omega0_de = 1. - (pba->Omega0_m + pba->Omega0_r + pba->Omega0_k);

  free(pvecback);
  free(pvecback_integration);

  return _SUCCESS_;

}

/**
 * Assign initial values to background integrated variables.
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input: pointer to background structure
 * @param pvecback             Input: vector of background quantities used as workspace
 * @param pvecback_integration Output: vector of background quantities to be integrated, returned with proper initial values
 * @return the error status
 */

int background_initial_conditions(
                                  struct precision *ppr,
                                  struct background *pba,
                                  double * pvecback, /* vector with argument pvecback[index_bg] (must be already allocated, normal format is sufficient) */
                                  double * pvecback_integration /* vector with argument pvecback_integration[index_bi] (must be already allocated with size pba->bi_size) */
                                  ) {

  /** Summary: */

  /** - define local variables */

  /* scale factor */
  double a;

  double rho_ncdm, p_ncdm, rho_ncdm_rel_tot=0.;
  double f,Omega_rad, rho_rad,T_ncdm,mu_ncdm, M, qmax;
  int counter,is_early_enough,n_ncdm,index_q;
  double scf_lambda;
  double rho_fld_today;
  double w_fld,dw_over_da_fld,integral_fld;

  /** - fix initial value of \f$ a \f$ */
  a = ppr->a_ini_over_a_today_default * pba->a_today;

  /**  If we have ncdm species, perhaps we need to start earlier
       than the standard value for the species to be relativistic.
       This could happen for some WDM models.
  */

  if (pba->has_ncdm == _TRUE_) {

    for (counter=0; counter < _MAX_IT_; counter++) {

      is_early_enough = _TRUE_;
      rho_ncdm_rel_tot = 0.;

      for (n_ncdm=0; n_ncdm<pba->N_ncdm; n_ncdm++) {

        if(pba->ncdm_background_distribution[n_ncdm] == _fermi_dirac_v2_ || pba->ncdm_background_distribution[n_ncdm] == _majoron_){
          class_call(interpolate_background_ncdm_distribution(pba,n_ncdm,pba->q_ncdm_bg[n_ncdm],pba->q_size_ncdm_bg[n_ncdm],1./a-1,pba->f_ncdm_bg[n_ncdm]),
          pba->error_message,
          pba->error_message);
          M = pba->m_ncdm_in_eV[n_ncdm];
          class_call(get_q_max(pba,n_ncdm,a,M,&qmax),
          pba->error_message,
          pba->error_message);
        }else{
          for(index_q = 0; index_q < pba->q_size_ncdm_bg[n_ncdm]; index_q++){
          pba->f_ncdm_bg[n_ncdm][index_q] = 1;//f_ncdm is already included in w_ncdm_bg in the case of standard neutrinos.
          }
          M = pba->M_ncdm[n_ncdm];
          qmax = 1;
        }



        class_call(background_ncdm_momenta(pba->q_ncdm_bg[n_ncdm],
                                           pba->w_ncdm_bg[n_ncdm],
                                           pba->f_ncdm_bg[n_ncdm],
                                           pba->q_size_ncdm_bg[n_ncdm],
                                           M,
                                           qmax,
                                           pba->factor_ncdm[n_ncdm],
                                           pba->a_today/a-1.0,
                                           n_ncdm,
                                           NULL,
                                           &rho_ncdm,
                                           &p_ncdm,
                                           NULL,
                                           NULL),
                   pba->error_message,
                   pba->error_message);
        rho_ncdm_rel_tot += 3.*p_ncdm;
        // printf("n_ncdm %d p_ncdm/rho_ncdm %e\n",n_ncdm,p_ncdm/rho_ncdm);
        if (fabs(p_ncdm/rho_ncdm-1./3.)>ppr->tol_ncdm_initial_w)
          is_early_enough = _FALSE_;
      }
      if (is_early_enough == _TRUE_)
        break;
      else
        a *= _SCALE_BACK_;
    }
    class_test(counter == _MAX_IT_,
               pba->error_message,
               "Search for initial scale factor a such that all ncdm species are relativistic failed.");
  }

  pvecback_integration[pba->index_bi_a] = a;

  /* Set initial values of {B} variables: */
  Omega_rad = pba->Omega0_g;
  if (pba->has_ur == _TRUE_)
    Omega_rad += pba->Omega0_ur;
  if (pba->has_idr == _TRUE_)
    Omega_rad += pba->Omega0_idr;
  rho_rad = Omega_rad*pow(pba->H0,2)/pow(a/pba->a_today,4);
  if (pba->has_ncdm == _TRUE_){
    /** - We must add the relativistic contribution from NCDM species */
    rho_rad += rho_ncdm_rel_tot;
  }
  if (pba->has_dcdm == _TRUE_){
    /* Remember that the critical density today in CLASS conventions is H0^2 */
    pvecback_integration[pba->index_bi_rho_dcdm] =
      pba->Omega_ini_dcdm*pba->H0*pba->H0*pow(pba->a_today/a,3);
    if (pba->background_verbose > 3)
      printf("Density is %g. a_today=%g. Omega_ini=%g\n",pvecback_integration[pba->index_bi_rho_dcdm],pba->a_today,pba->Omega_ini_dcdm);
  }

  if (pba->has_dr == _TRUE_){
    if (pba->has_dcdm == _TRUE_){
      /**  - f is the critical density fraction of DR. The exact solution is:
       *
       * `f = -Omega_rad+pow(pow(Omega_rad,3./2.)+0.5*pow(a/pba->a_today,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3),2./3.);`
       *
       * but it is not numerically stable for very small f which is always the case.
       * Instead we use the Taylor expansion of this equation, which is equivalent to
       * ignoring f(a) in the Hubble rate.
       */
      f = 1./3.*pow(a/pba->a_today,6)*pvecback_integration[pba->index_bi_rho_dcdm]*pba->Gamma_dcdm/pow(pba->H0,3)/sqrt(Omega_rad);
      pvecback_integration[pba->index_bi_rho_dr] = f*pba->H0*pba->H0/pow(a/pba->a_today,4);
    }
    else{
      /** There is also a space reserved for a future case where dr is not sourced by dcdm */
      pvecback_integration[pba->index_bi_rho_dr] = 0.0;
    }
  }

  if (pba->has_fld == _TRUE_){

    /* rho_fld today */
    rho_fld_today = pba->Omega0_fld * pow(pba->H0,2);

    /* integrate rho_fld(a) from a_ini to a_0, to get rho_fld(a_ini) given rho_fld(a0) */
    class_call(background_w_fld(pba,a,&w_fld,&dw_over_da_fld,&integral_fld), pba->error_message, pba->error_message);

    /* Note: for complicated w_fld(a) functions with no simple
       analytic integral, this is the place were you should compute
       numerically the simple 1d integral [int_{a_ini}^{a_0} 3
       [(1+w_fld)/a] da] (e.g. with the Romberg method?) instead of
       calling background_w_fld */

    /* rho_fld at initial time */
    pvecback_integration[pba->index_bi_rho_fld] = rho_fld_today * exp(integral_fld);

  }

  /** - Fix initial value of \f$ \phi, \phi' \f$
   * set directly in the radiation attractor => fixes the units in terms of rho_ur
   *
   * TODO:
   * - There seems to be some small oscillation when it starts.
   * - Check equations and signs. Sign of phi_prime?
   * - is rho_ur all there is early on?
   */
  if(pba->has_scf == _TRUE_){
    scf_lambda = pba->scf_parameters[0];
    if(pba->attractor_ic_scf == _TRUE_){
      pvecback_integration[pba->index_bi_phi_scf] = -1/scf_lambda*
        log(rho_rad*4./(3*pow(scf_lambda,2)-12))*pba->phi_ini_scf;
      if (3.*pow(scf_lambda,2)-12. < 0){
        /** - --> If there is no attractor solution for scf_lambda, assign some value. Otherwise would give a nan.*/
    	pvecback_integration[pba->index_bi_phi_scf] = 1./scf_lambda;//seems to the work
        if (pba->background_verbose > 0)
          printf(" No attractor IC for lambda = %.3e ! \n ",scf_lambda);
      }
      pvecback_integration[pba->index_bi_phi_prime_scf] = 2*pvecback_integration[pba->index_bi_a]*
        sqrt(V_scf(pba,pvecback_integration[pba->index_bi_phi_scf]))*pba->phi_prime_ini_scf;
    }
    else{
      printf("Not using attractor initial conditions\n");
      /** - --> If no attractor initial conditions are assigned, gets the provided ones. */
      pvecback_integration[pba->index_bi_phi_scf] = pba->phi_ini_scf;
      pvecback_integration[pba->index_bi_phi_prime_scf] = pba->phi_prime_ini_scf;
    }
    class_test(!isfinite(pvecback_integration[pba->index_bi_phi_scf]) ||
               !isfinite(pvecback_integration[pba->index_bi_phi_scf]),
               pba->error_message,
               "initial phi = %e phi_prime = %e -> check initial conditions",
               pvecback_integration[pba->index_bi_phi_scf],
               pvecback_integration[pba->index_bi_phi_scf]);
  }

  /* Infer pvecback from pvecback_integration */


  class_call(background_functions(pba, pvecback_integration, pba->normal_info, pvecback),
             pba->error_message,
             pba->error_message);

  /* Just checking that our initial time indeed is deep enough in the radiation
     dominated regime */
  class_test(fabs(pvecback[pba->index_bg_Omega_r]-1.) > ppr->tol_initial_Omega_r,
             pba->error_message,
             "Omega_r = %e, not close enough to 1. Decrease a_ini_over_a_today_default in order to start from radiation domination.",
             pvecback[pba->index_bg_Omega_r]);

  /** - compute initial proper time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ t=1/(2H) \f$ (good
      approximation for most purposes) */

  class_test(pvecback[pba->index_bg_H] <= 0.,
             pba->error_message,
             "H = %e instead of strictly positive",pvecback[pba->index_bg_H]);

  pvecback_integration[pba->index_bi_time] = 1./(2.* pvecback[pba->index_bg_H]);

  /** - compute initial conformal time, assuming radiation-dominated
      universe since Big Bang and therefore \f$ \tau=1/(aH) \f$
      (good approximation for most purposes) */
  pvecback_integration[pba->index_bi_tau] = 1./(a * pvecback[pba->index_bg_H]);

  /** - compute initial sound horizon, assuming \f$ c_s=1/\sqrt{3} \f$ initially */
  pvecback_integration[pba->index_bi_rs] = pvecback_integration[pba->index_bi_tau]/sqrt(3.);

  /** - set initial value of D and D' in RD. D will be renormalised later, but D' must be correct. */
  pvecback_integration[pba->index_bi_D] = a;
  pvecback_integration[pba->index_bi_D_prime] = 2*pvecback_integration[pba->index_bi_D]*pvecback[pba->index_bg_H];

  return _SUCCESS_;

}

/**
 * Find the time of radiation/matter equality and store characteristic
 * quantitites at that time in the background structure..
 *
 * @param ppr                  Input: pointer to precision structure
 * @param pba                  Input/Output: pointer to background structure
 * @return the error status
 */

int background_find_equality(
                             struct precision *ppr,
                             struct background *pba) {

  double Omega_m_over_Omega_r=0.;
  int index_tau_minus = 0;
  int index_tau_plus = pba->bt_size-1;
  int index_tau_mid = 0;
  double tau_minus,tau_plus,tau_mid=0.;
  double * pvecback;

  /* first bracket the right tau value between two consecutive indices in the table */

  while ((index_tau_plus - index_tau_minus) > 1) {

    index_tau_mid = (int)(0.5*(index_tau_plus+index_tau_minus));

    Omega_m_over_Omega_r = pba->background_table[index_tau_mid*pba->bg_size+pba->index_bg_Omega_m]
      /pba->background_table[index_tau_mid*pba->bg_size+pba->index_bg_Omega_r];

    if (Omega_m_over_Omega_r > 1)
      index_tau_plus = index_tau_mid;
    else
      index_tau_minus = index_tau_mid;

  }

  /* then get a better estimate within this range */

  tau_minus = pba->tau_table[index_tau_minus];
  tau_plus =  pba->tau_table[index_tau_plus];

  class_alloc(pvecback,pba->bg_size*sizeof(double),pba->error_message);

  while ((tau_plus - tau_minus) > ppr->tol_tau_eq) {

    tau_mid = 0.5*(tau_plus+tau_minus);

    class_call(background_at_tau(pba,tau_mid,pba->long_info,pba->inter_closeby,&index_tau_minus,pvecback),
               pba->error_message,
               pba->error_message);

    Omega_m_over_Omega_r = pvecback[pba->index_bg_Omega_m]/pvecback[pba->index_bg_Omega_r];

    if (Omega_m_over_Omega_r > 1)
      tau_plus = tau_mid;
    else
      tau_minus = tau_mid;

  }

  pba->a_eq = pvecback[pba->index_bg_a];
  pba->H_eq = pvecback[pba->index_bg_H];
  pba->z_eq = pba->a_today/pba->a_eq -1.;
  pba->tau_eq = tau_mid;

  if (pba->background_verbose > 0) {
    printf(" -> radiation/matter equality at z = %f\n",pba->z_eq);
    printf("    corresponding to conformal time = %f Mpc\n",pba->tau_eq);
  }

  free(pvecback);

  return _SUCCESS_;

}


/**
 * Subroutine for formatting background output
 *
 */

int background_output_titles(struct background * pba,
                             char titles[_MAXTITLESTRINGLENGTH_]
                             ){

  /** - Length of the column title should be less than _OUTPUTPRECISION_+6
      to be indented correctly, but it can be as long as . */
  int n;
  char tmp[40];

  class_store_columntitle(titles,"z",_TRUE_);
  class_store_columntitle(titles,"proper time [Gyr]",_TRUE_);
  class_store_columntitle(titles,"conf. time [Mpc]",_TRUE_);
  class_store_columntitle(titles,"H [1/Mpc]",_TRUE_);
  class_store_columntitle(titles,"comov. dist.",_TRUE_);
  class_store_columntitle(titles,"ang.diam.dist.",_TRUE_);
  class_store_columntitle(titles,"lum. dist.",_TRUE_);
  class_store_columntitle(titles,"comov.snd.hrz.",_TRUE_);
  class_store_columntitle(titles,"(.)rho_g",_TRUE_);
  class_store_columntitle(titles,"(.)rho_b",_TRUE_);
  class_store_columntitle(titles,"(.)rho_cdm",pba->has_cdm);
  if (pba->has_ncdm == _TRUE_){
    for (n=0; n<pba->N_ncdm; n++){
      sprintf(tmp,"(.)rho_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      sprintf(tmp,"(.)p_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      sprintf(tmp,"(.)pseudo_p_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      sprintf(tmp,"(.)T_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      sprintf(tmp,"(.)dT_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      sprintf(tmp,"(.)mu_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
      sprintf(tmp,"(.)dmu_ncdm[%d]",n);
      class_store_columntitle(titles,tmp,_TRUE_);
    }
  }
  class_store_columntitle(titles,"(.)rho_lambda",pba->has_lambda);
  class_store_columntitle(titles,"(.)rho_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)w_fld",pba->has_fld);
  class_store_columntitle(titles,"(.)rho_ur",pba->has_ur);
  class_store_columntitle(titles,"(.)rho_idr",pba->has_idr);
  class_store_columntitle(titles,"(.)rho_idm_dr",pba->has_idm_dr);
  class_store_columntitle(titles,"(.)rho_crit",_TRUE_);
  class_store_columntitle(titles,"(.)rho_dcdm",pba->has_dcdm);
  class_store_columntitle(titles,"(.)rho_dr",pba->has_dr);

  class_store_columntitle(titles,"(.)rho_scf",pba->has_scf);
  class_store_columntitle(titles,"(.)p_scf",pba->has_scf);
  class_store_columntitle(titles,"(.)p_prime_scf",pba->has_scf);
  class_store_columntitle(titles,"phi_scf",pba->has_scf);
  class_store_columntitle(titles,"phi'_scf",pba->has_scf);
  class_store_columntitle(titles,"V_scf",pba->has_scf);
  class_store_columntitle(titles,"V'_scf",pba->has_scf);
  class_store_columntitle(titles,"V''_scf",pba->has_scf);

  class_store_columntitle(titles,"(.)rho_tot",_TRUE_);
  class_store_columntitle(titles,"(.)p_tot",_TRUE_);
  class_store_columntitle(titles,"(.)p_tot_prime",_TRUE_);

  class_store_columntitle(titles,"gr.fac. D",_TRUE_);
  class_store_columntitle(titles,"gr.fac. f",_TRUE_);

  return _SUCCESS_;
}

int background_output_data(
                           struct background *pba,
                           int number_of_titles,
                           double *data){
  int index_tau, storeidx, n;
  double *dataptr, *pvecback;

  /** Stores quantities */
  for (index_tau=0; index_tau<pba->bt_size; index_tau++){
    dataptr = data + index_tau*number_of_titles;
    pvecback = pba->background_table + index_tau*pba->bg_size;
    storeidx = 0;

    class_store_double(dataptr,pba->a_today/pvecback[pba->index_bg_a]-1.,_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_time]/_Gyr_over_Mpc_,_TRUE_,storeidx);
    class_store_double(dataptr,pba->conformal_age-pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_H],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_conf_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_ang_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_lum_distance],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rs],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_g],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_b],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_cdm],pba->has_cdm,storeidx);
    if (pba->has_ncdm == _TRUE_){
      for (n=0; n<pba->N_ncdm; n++){
        class_store_double(dataptr,pvecback[pba->index_bg_rho_ncdm1+n],_TRUE_,storeidx);
        class_store_double(dataptr,pvecback[pba->index_bg_p_ncdm1+n],_TRUE_,storeidx);
        class_store_double(dataptr,pvecback[pba->index_bg_pseudo_p_ncdm1+n],_TRUE_,storeidx);
        class_store_double(dataptr,pvecback[pba->index_bg_T_ncdm1+n],_TRUE_,storeidx);
        class_store_double(dataptr,pvecback[pba->index_bg_dT_ncdm1+n],_TRUE_,storeidx);
        class_store_double(dataptr,pvecback[pba->index_bg_Mu_ncdm1+n],_TRUE_,storeidx);
        class_store_double(dataptr,pvecback[pba->index_bg_dMu_ncdm1+n],_TRUE_,storeidx);
      }
    }
    class_store_double(dataptr,pvecback[pba->index_bg_rho_lambda],pba->has_lambda,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_fld],pba->has_fld,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_w_fld],pba->has_fld,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_ur],pba->has_ur,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_idr],pba->has_idr,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_idm_dr],pba->has_idm_dr,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_crit],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dcdm],pba->has_dcdm,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_rho_dr],pba->has_dr,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_rho_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_prime_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_phi_prime_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_V_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_dV_scf],pba->has_scf,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_ddV_scf],pba->has_scf,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_rho_tot],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_tot],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_p_tot_prime],_TRUE_,storeidx);

    class_store_double(dataptr,pvecback[pba->index_bg_D],_TRUE_,storeidx);
    class_store_double(dataptr,pvecback[pba->index_bg_f],_TRUE_,storeidx);
  }

  return _SUCCESS_;
}


/**
 * Subroutine evaluating the derivative with respect to conformal time
 * of quantities which are integrated (a, t, etc).
 *
 * This is one of the few functions in the code which is passed to
 * the generic_integrator() routine.  Since generic_integrator()
 * should work with functions passed from various modules, the format
 * of the arguments is a bit special:
 *
 * - fixed input parameters and workspaces are passed through a generic
 * pointer. Here, this is just a pointer to the background structure
 * and to a background vector, but generic_integrator() doesn't know
 * its fine structure.
 *
 * - the error management is a bit special: errors are not written as
 * usual to pba->error_message, but to a generic error_message passed
 * in the list of arguments.
 *
 * @param tau                      Input: conformal time
 * @param y                        Input: vector of variable
 * @param dy                       Output: its derivative (already allocated)
 * @param parameters_and_workspace Input: pointer to fixed parameters (e.g. indices)
 * @param error_message            Output: error message
 */
int background_derivs(
                      double tau,
                      double* y, /* vector with argument y[index_bi] (must be already allocated with size pba->bi_size) */
                      double* dy, /* vector with argument dy[index_bi]
                                     (must be already allocated with
                                     size pba->bi_size) */
                      void * parameters_and_workspace,
                      ErrorMsg error_message
                      ) {

  /** Summary: */

  /** - define local variables */

  struct background_parameters_and_workspace * pbpaw;
  struct background * pba;
  double * pvecback, a, H, rho_M;

  pbpaw = parameters_and_workspace;
  pba =  pbpaw->pba;
  pvecback = pbpaw->pvecback;

  /** - calculate functions of \f$ a \f$ with background_functions() */
  class_call(background_functions(pba, y, pba->normal_info, pvecback),
             pba->error_message,
             error_message);

  /** - Short hand notation */
  a = y[pba->index_bi_a];
  H = pvecback[pba->index_bg_H];

  /** - calculate \f$ a'=a^2 H \f$ */
  dy[pba->index_bi_a] = y[pba->index_bi_a] * y[pba->index_bi_a] * pvecback[pba->index_bg_H];

  /** - calculate \f$ t' = a \f$ */
  dy[pba->index_bi_time] = y[pba->index_bi_a];

  class_test(pvecback[pba->index_bg_rho_g] <= 0.,
             error_message,
             "rho_g = %e instead of strictly positive",pvecback[pba->index_bg_rho_g]);

  /** - calculate \f$ rs' = c_s \f$*/
  dy[pba->index_bi_rs] = 1./sqrt(3.*(1.+3.*pvecback[pba->index_bg_rho_b]/4./pvecback[pba->index_bg_rho_g]))*sqrt(1.-pba->K*y[pba->index_bi_rs]*y[pba->index_bi_rs]); // TBC: curvature correction

  /** - solve second order growth equation  \f$ [D''(\tau)=-aHD'(\tau)+3/2 a^2 \rho_M D(\tau) \f$ */
  rho_M = pvecback[pba->index_bg_rho_b];
  if (pba->has_cdm)
    rho_M += pvecback[pba->index_bg_rho_cdm];
  if (pba->has_idm_dr)
    rho_M += pvecback[pba->index_bg_rho_idm_dr];

  dy[pba->index_bi_D] = y[pba->index_bi_D_prime];
  dy[pba->index_bi_D_prime] = -a*H*y[pba->index_bi_D_prime] + 1.5*a*a*rho_M*y[pba->index_bi_D];

  if (pba->has_dcdm == _TRUE_){
    /** - compute dcdm density \f$ \rho' = -3aH \rho - a \Gamma \rho \f$*/
    dy[pba->index_bi_rho_dcdm] = -3.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*y[pba->index_bi_rho_dcdm]-
      y[pba->index_bi_a]*pba->Gamma_dcdm*y[pba->index_bi_rho_dcdm];
  }

  if ((pba->has_dcdm == _TRUE_) && (pba->has_dr == _TRUE_)){
    /** - Compute dr density \f$ \rho' = -4aH \rho - a \Gamma \rho \f$ */
    dy[pba->index_bi_rho_dr] = -4.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*y[pba->index_bi_rho_dr]+
      y[pba->index_bi_a]*pba->Gamma_dcdm*y[pba->index_bi_rho_dcdm];
  }

  if (pba->has_fld == _TRUE_) {
    /** - Compute fld density \f$ \rho' = -3aH (1+w_{fld}(a)) \rho \f$ */
    dy[pba->index_bi_rho_fld] = -3.*y[pba->index_bi_a]*pvecback[pba->index_bg_H]*(1.+pvecback[pba->index_bg_w_fld])*y[pba->index_bi_rho_fld];
  }

  if (pba->has_scf == _TRUE_){
    /** - Scalar field equation: \f$ \phi'' + 2 a H \phi' + a^2 dV = 0 \f$  (note H is wrt cosmic time) */
    dy[pba->index_bi_phi_scf] = y[pba->index_bi_phi_prime_scf];
    dy[pba->index_bi_phi_prime_scf] = - y[pba->index_bi_a]*
      (2*pvecback[pba->index_bg_H]*y[pba->index_bi_phi_prime_scf]
       + y[pba->index_bi_a]*dV_scf(pba,y[pba->index_bi_phi_scf])) ;
  }

  return _SUCCESS_;

}

/**
 * Scalar field potential and its derivatives with respect to the field _scf
 * For Albrecht & Skordis model: 9908085
 * - \f$ V = V_{p_{scf}}*V_{e_{scf}} \f$
 * - \f$ V_e =  \exp(-\lambda \phi) \f$ (exponential)
 * - \f$ V_p = (\phi - B)^\alpha + A \f$ (polynomial bump)
 *
 * TODO:
 * - Add some functionality to include different models/potentials (tuning would be difficult, though)
 * - Generalize to Kessence/Horndeski/PPF and/or couplings
 * - A default module to numerically compute the derivatives when no analytic functions are given should be added.
 * - Numerical derivatives may further serve as a consistency check.
 *
 */

/**
 *
 * The units of phi, tau in the derivatives and the potential V are the following:
 * - phi is given in units of the reduced Planck mass \f$ m_{pl} = (8 \pi G)^{(-1/2)}\f$
 * - tau in the derivative is given in units of Mpc.
 * - the potential \f$ V(\phi) \f$ is given in units of \f$ m_{pl}^2/Mpc^2 \f$.
 * With this convention, we have
 * \f$ \rho^{class} = (8 \pi G)/3 \rho^{physical} = 1/(3 m_{pl}^2) \rho^{physical} = 1/3 * [ 1/(2a^2) (\phi')^2 + V(\phi) ] \f$
 and \f$ \rho^{class} \f$ has the proper dimension \f$ Mpc^-2 \f$.
*/

double V_e_scf(struct background *pba,
               double phi
               ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return  exp(-scf_lambda*phi);
}

double dV_e_scf(struct background *pba,
                double phi
                ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return -scf_lambda*V_scf(pba,phi);
}

double ddV_e_scf(struct background *pba,
                 double phi
                 ) {
  double scf_lambda = pba->scf_parameters[0];
  //  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  //  double scf_B      = pba->scf_parameters[3];

  return pow(-scf_lambda,2)*V_scf(pba,phi);
}


/** parameters and functions for the polynomial coefficient
 * \f$ V_p = (\phi - B)^\alpha + A \f$(polynomial bump)
 *
 * double scf_alpha = 2;
 *
 * double scf_B = 34.8;
 *
 * double scf_A = 0.01; (values for their Figure 2)
 */

double V_p_scf(
               struct background *pba,
               double phi) {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  pow(phi - scf_B,  scf_alpha) +  scf_A;
}

double dV_p_scf(
                struct background *pba,
                double phi) {

  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return   scf_alpha*pow(phi -  scf_B,  scf_alpha - 1);
}

double ddV_p_scf(
                 struct background *pba,
                 double phi) {
  //  double scf_lambda = pba->scf_parameters[0];
  double scf_alpha  = pba->scf_parameters[1];
  //  double scf_A      = pba->scf_parameters[2];
  double scf_B      = pba->scf_parameters[3];

  return  scf_alpha*(scf_alpha - 1.)*pow(phi -  scf_B,  scf_alpha - 2);
}

/** Fianlly we can obtain the overall potential \f$ V = V_p*V_e \f$
 */

double V_scf(
             struct background *pba,
             double phi) {
  return  V_e_scf(pba,phi)*V_p_scf(pba,phi);
}

double dV_scf(
              struct background *pba,
              double phi) {
  return dV_e_scf(pba,phi)*V_p_scf(pba,phi) + V_e_scf(pba,phi)*dV_p_scf(pba,phi);
}

double ddV_scf(
               struct background *pba,
               double phi) {
  return ddV_e_scf(pba,phi)*V_p_scf(pba,phi) + 2*dV_e_scf(pba,phi)*dV_p_scf(pba,phi) + V_e_scf(pba,phi)*ddV_p_scf(pba,phi);
}

/**
 * Function outputting the fractions Omega of the total critical density
 * today, and also the reduced fractions omega=Omega*h*h
 *
 * It also prints the total budgets of non-relativistic, relativistic,
 * and other contents, and of the total
 *
 * @param pba                      Input: Pointer to background structure
 * @return the error status
 */

int background_output_budget(
                             struct background* pba
                             ) {

  double budget_matter, budget_radiation, budget_other,budget_neutrino;
  int index_ncdm;

  budget_matter = 0;
  budget_radiation = 0;
  budget_other = 0;
  budget_neutrino = 0;

  //The name for the _class_print_species_ macro can be at most 30 characters total
  if(pba->background_verbose > 1){

    printf(" ---------------------------- Budget equation ----------------------- \n");

    printf(" ---> Nonrelativistic Species \n");
    _class_print_species_("Bayrons",b);
    budget_matter+=pba->Omega0_b;
    if(pba->has_cdm){
      _class_print_species_("Cold Dark Matter",cdm);
      budget_matter+=pba->Omega0_cdm;
    }
    if(pba->has_idm_dr){
      _class_print_species_("Interacting Dark Matter - DR ",idm_dr);
      budget_matter+=pba->Omega0_idm_dr;
    }
    if(pba->has_dcdm){
      _class_print_species_("Decaying Cold Dark Matter",dcdm);
      budget_matter+=pba->Omega0_dcdm;
    }


    printf(" ---> Relativistic Species \n");
    _class_print_species_("Photons",g);
    budget_radiation+=pba->Omega0_g;
    if(pba->has_ur){
      _class_print_species_("Ultra-relativistic relics",ur);
      budget_radiation+=pba->Omega0_ur;
    }
    if(pba->has_dr){
      _class_print_species_("Dark Radiation (from decay)",dr);
      budget_radiation+=pba->Omega0_dr;
    }
    if(pba->has_idr){
      _class_print_species_("Interacting Dark Radiation",idr);
      budget_radiation+=pba->Omega0_idr;
    }

    if(pba->N_ncdm > 0){
      printf(" ---> Massive Neutrino Species \n");
    }
    if(pba->N_ncdm > 0){
      for(index_ncdm=0;index_ncdm<pba->N_ncdm;++index_ncdm){
        printf("-> %-26s%-4d Omega = %-15g , omega = %-15g\n","Neutrino Species Nr.",index_ncdm+1,pba->Omega0_ncdm[index_ncdm],pba->Omega0_ncdm[index_ncdm]*pba->h*pba->h);
        budget_neutrino+=pba->Omega0_ncdm[index_ncdm];
      }
    }

    if(pba->has_lambda || pba->has_fld || pba->has_scf || pba->has_curvature){
      printf(" ---> Other Content \n");
    }
    if(pba->has_lambda){
      _class_print_species_("Cosmological Constant",lambda);
      budget_other+=pba->Omega0_lambda;
    }
    if(pba->has_fld){
      _class_print_species_("Dark Energy Fluid",fld);
      budget_other+=pba->Omega0_fld;
    }
    if(pba->has_scf){
      _class_print_species_("Scalar Field",scf);
      budget_other+=pba->Omega0_scf;
    }
    if(pba->has_curvature){
      _class_print_species_("Spatial Curvature",k);
      budget_other+=pba->Omega0_k;
    }

    printf(" ---> Total budgets \n");
    printf(" Radiation                        Omega = %-15g , omega = %-15g \n",budget_radiation,budget_radiation*pba->h*pba->h);
    printf(" Non-relativistic                 Omega = %-15g , omega = %-15g \n",budget_matter,budget_matter*pba->h*pba->h);
    if(pba->N_ncdm > 0){
      printf(" Neutrinos                        Omega = %-15g , omega = %-15g \n",budget_neutrino,budget_neutrino*pba->h*pba->h);
    }
    if(pba->has_lambda || pba->has_fld || pba->has_scf || pba->has_curvature){
      printf(" Other Content                    Omega = %-15g , omega = %-15g \n",budget_other,budget_other*pba->h*pba->h);
    }
    printf(" TOTAL                            Omega = %-15g , omega = %-15g \n",budget_radiation+budget_matter+budget_neutrino+budget_other,(budget_radiation+budget_matter+budget_neutrino+budget_other)*pba->h*pba->h);

    printf(" -------------------------------------------------------------------- \n");
  }

  return _SUCCESS_;
}


/**
 * SJW: New function
 */

int background_MB_approx(
//                            /* Only calculate for non-NULL pointers: */
                            struct background *pba,
                            int *lenIndx) {
  double Mnu, mMaj;
  double GammaEff, GammaPhi;
  double tcur, presNu, presMaj;
  double rhoMaj, nMaj, rhoNu, nNu;
  double muMh, muMstore;
  double Hub, h_cmb, H_preF, h_mat;
  double degNu=6, degMaj=1, neff;
  double k1[5], k2[5], k3[5], k4[5], tmajH, tnuH, muNh;
  int sigP=100, indx, worked;
  double sigSt, ehold1, ehold0, trap0, trap1;
  double maxBndMaj, maxBndNu, ThreshJump;
  double zstart, zend, zhold, delT;
  bool lower_sve_indx=false, shrinkDT=false, stop_loop=false, linearMu=false, linearMuN=false, linearT=false;
  int numT=500000, sve_indx=100;
  *lenIndx=0;

  GammaPhi = pba->Gamma_phi[0]; // Not yet generalized to deal with multiple neutrinos....
  for (int n_ncdm=0; n_ncdm < pba->N_ncdm; n_ncdm++) {
        // if(pba->ncdm_background_distribution[n_ncdm]==_majoron_){mMaj = pba->m_ncdm_in_eV[n_ncdm];}
        if(n_ncdm == pba->entry_is_M_phi){mMaj = pba->m_ncdm_in_eV[n_ncdm];}
        else{Mnu = pba->m_ncdm_in_eV[n_ncdm];}
  }

  GammaEff = pba->Gamma_eff_ncdmphi[0];
  if (GammaEff > 1e1) { GammaPhi = GammaPhi * 10. / GammaEff;}

  ThreshJump = mMaj * 2.;
  zstart = 100. * mMaj / (pba->T_cmb * 8.6e-5) *  1.39578 - 1;
  zend = mMaj / (pba->T_cmb * 8.6e-5) - 1.;
  delT = (zend - zstart) / numT * 10;
  double muNstore, TNstore, TMstore, TcmbStore, zstore;
  zstore = zstart;
  TcmbStore = pba->T_cmb * (1. + zstart)  * 8.6e-5; // eV
  TNstore = pba->T_cmb * (1. + zstart) / 1.39578 * 8.6e-5; // eV

  TMstore = 1e-2 * TNstore; // eV
  muNstore = -1e-3 * TNstore; // eV
  muMstore = - 10. * TMstore; // eV
  H_preF = 8.*_PI_/3. * 6.707e-57; // eV^-2
  indx = 0;

  while (!stop_loop){
    indx+=1;

    if ((TcmbStore < ThreshJump)&&(!shrinkDT)){
    delT /= 10.;
    if (ThreshJump == mMaj){
        ThreshJump = 0.6 * mMaj;
        }
    else {shrinkDT=true;}
    }


    // evaluated at z_n -- get k_1
    tmajH = TMstore;
    tnuH = TNstore;
    muMh = muMstore;
    muNh = muNstore;
    tcur = TcmbStore;
    zhold = zstore;

    worked = RK_Eval(pba, GammaPhi, zhold, tmajH, tnuH, muMh, muNh, mMaj, tcur, Mnu, k1);
    if (linearMu) {k1[2] *= muMh;}
    if (linearMuN) {k1[3] *= muNh;}
    //if (linearT) {k1[0] *= tmajH;}
    //printf("K1: %e \t %e \t %e \t %e \t %e \n", k1[0], k1[1], k1[2], k1[3], k1[4]);

    // get k_2, f(z + h/2, y + h k_1/2)
    //tmajH = exp(log(TMstore) + delT * k1[0] / 2);
    if (!linearT) {tmajH = exp(log(TMstore) + delT * k1[0] / 2);}
    else {tmajH = TMstore +  delT * k1[0] / 2;}
    tnuH = exp(log(TNstore)   + delT * k1[1] / 2);
    if (!linearMu) {muMh = -exp(log(-muMstore) + delT * k1[2] / 2);}
    else {muMh = muMstore +  delT * k1[2] / 2;}
    if (!linearMuN) {muNh = -exp(log(-muNstore)  + delT * k1[3] / 2);}
    else {muNh = muNstore +  delT * k1[3] / 2;}

    tcur = TcmbStore + delT * k1[4] / 2;
    zhold = zstore + delT / 2.;

    worked = RK_Eval(pba, GammaPhi, zhold, tmajH, tnuH, muMh, muNh, mMaj, tcur, Mnu, k2);
    if (linearMu) {k2[2] *= muMh;}
    if (linearMuN) {k2[3] *= muNh;}
    //if (linearT) {k2[0] *= tmajH;}
    //printf("K2: %e \t %e \t %e \t %e \t %e \n", k2[0], k2[1], k2[2], k2[3], k2[4]);

    // get k_3, f(z + h/2, y + h k_2/2)
    //tmajH = exp(log(TMstore) + delT * k2[0] / 2);
    if (!linearT) {tmajH = exp(log(TMstore) + delT * k2[0] / 2);}
    else {tmajH = TMstore +  delT * k2[0] / 2;}
    tnuH = exp(log(TNstore)   + delT * k2[1] / 2);
    if (!linearMu) {muMh = -exp(log(-muMstore) + delT * k2[2] / 2);}
    else {muMh = muMstore +  delT * k2[2] / 2;}
    if (!linearMuN) {muNh = -exp(log(-muNstore)  + delT * k2[3] / 2);}
    else {muNh = muNstore +  delT * k2[3] / 2;}

    tcur = TcmbStore + delT * k2[4] / 2;
    zhold = zstore + delT / 2.;

    worked = RK_Eval(pba, GammaPhi, zhold, tmajH, tnuH, muMh, muNh, mMaj, tcur, Mnu, k3);
    if (linearMu) {k3[2] *= muMh;}
    if (linearMuN) {k3[3] *= muNh;}
    //if (linearT) {k3[0] *= tmajH;}
    //printf("K3: %e \t %e \t %e \t %e \t %e \n", k3[0], k3[1], k3[2], k3[3], k3[4]);

    // get k_4, f(z + h, y + h k_3)
    //tmajH = exp(log(TMstore) + delT * k3[0]);
    if (!linearT) {tmajH = exp(log(TMstore) + delT * k3[0]);}
    else {tmajH = TMstore +  delT * k3[0];}
    tnuH = exp(log(TNstore)   + delT * k3[1]);
    if (!linearMu) {muMh = -exp(log(-muMstore) + delT * k3[2]);}
    else {muMh = muMstore +  delT * k3[2];}
    if (!linearMu) {muNh = -exp(log(-muNstore)  + delT * k3[3]);}
    else {muNh = muNstore +  delT * k3[3];}

    tcur = TcmbStore  + delT * k3[4];
    zhold = zstore + delT;
    worked = RK_Eval(pba, GammaPhi, zhold, tmajH, tnuH, muMh, muNh, mMaj, tcur, Mnu, k4);
    if (linearMu) {k4[2] *= muMh;}
    if (linearMuN) {k4[3] *= muNh;}
    //if (linearT) {k4[0] *= tmajH;}
    //printf("K4: %e \t %e \t %e \t %e \t %e \n", k4[0], k4[1], k4[2], k4[3], k4[4]);


    // Fill new variables
    if (!linearT) {TMstore = exp(log(TMstore) + delT / 6. * (k1[0] + 2. * (k2[0] + k3[0]) + k4[0]));}
    else {TMstore = TMstore +  delT / 6. * (k1[0] + 2. * (k2[0] + k3[0]) + k4[0]);}
    if (!linearMu) {muMstore = - exp(log(-muMstore) +  delT / 6. * (k1[2] + 2. * (k2[2] + k3[2]) + k4[2]));}
    else {muMstore = muMstore +  delT / 6. * (k1[2] + 2. * (k2[2] + k3[2]) + k4[2]);}
    if (!linearMuN) {muNstore = -exp(log(-muNstore) + delT / 6. * (k1[3] + 2. * (k2[3] + k3[3]) + k4[3]));}
    else {muNstore = muNstore +  delT / 6. * (k1[3] + 2. * (k2[3] + k3[3]) + k4[3]);}

    TNstore = exp(log(TNstore) + delT / 6. * (k1[1] + 2. * (k2[1] + k3[1]) + k4[1]));
    TcmbStore = TcmbStore + delT / 6. * (k1[4] + 2. * (k2[4] + k3[4]) + k4[4]);
    zstore = zstore + delT;

    if ((!linearMu) && (muMstore > -5e-3) ){
    //printf("Switching to linear.... \t %e \n", muMstore);
    linearMu = true;
//    linearT = true;
    }
    if ((!linearMuN) && (fabs(muNstore) < 1e-4) ){linearMuN = true;    }


    if (!isfinite(TMstore)||!isfinite(TNstore)||!isfinite(muNstore)){
        //printf("Failure.... Nans or infs appearing... \n");
        //printf("%e \t %e \t %e \n", TMstore, TNstore, muNstore);
        break;
    }


    if (!lower_sve_indx&&(TcmbStore < 0.8 * mMaj)) {
    sve_indx = 5;
    lower_sve_indx = true;
    }
    // How often commpute energy densities, pressure, etc.
    if (indx % sve_indx == 0){
        tmajH = TMstore;
        tnuH = TNstore;
        muMh = muMstore;
        muNh = muNstore;
        tcur = TcmbStore;
        zhold = zstore;

        h_cmb = pow(_PI_, 2) / 15. * pow(tcur, 4); // energy density in ev^4
        h_mat = (pba->Omega0_b + pba->Omega0_cdm) * pow((1. + zhold), 3) * 8.0835e-11 * pow(2.998e10 * 6.58e-16, 3) * pow(pba->h, 2); // eV^4

        maxBndMaj = 20. * tmajH;
        if (maxBndMaj < (3. * mMaj)){maxBndMaj=3.*mMaj;}
        nMaj = 0.;
        rhoMaj = 0.;
        presMaj=0.;
        // Maj integration
        sigSt = (double)((maxBndMaj - mMaj) / sigP);
        for (int ii=0; ii < (sigP-1); ii++) {
            ehold0 = mMaj + sigSt * (double)ii;
            ehold1 = mMaj + sigSt * (double)(ii+1);
            if (ii==0) {ehold0 += 1e-20;}
            // n
            trap0 = ehold0 * sqrt(ehold0*ehold0 - mMaj*mMaj) /  ( exp((ehold0 - muMh) / tmajH) - 1.) ;
            trap1 = ehold1 * sqrt(ehold1*ehold1 - mMaj*mMaj) /  ( exp((ehold1 - muMh) / tmajH) - 1.);
            nMaj += 0.5 * sigSt * (trap0 + trap1)  / (2 * _PI_ * _PI_);
            // rho
            trap0 = ehold0 * ehold0 * sqrt(ehold0*ehold0 - mMaj*mMaj) /  ( exp((ehold0 - muMh) / tmajH) - 1.);
            trap1 = ehold1 * ehold1 * sqrt(ehold1*ehold1 - mMaj*mMaj) /  ( exp((ehold1 - muMh) / tmajH) - 1.) ;
            rhoMaj += 0.5 * sigSt * (trap0 + trap1)  / (2 * _PI_ * _PI_);
            // Pres
            trap0 =   pow(ehold0*ehold0 - mMaj*mMaj, 1.5) /  ( exp((ehold0 - muMh) / tmajH) - 1.) ;
            trap1 =  pow(ehold1*ehold1 - mMaj*mMaj, 1.5) /  ( exp((ehold1 - muMh) / tmajH) - 1.) ;
            presMaj += degMaj * 0.5 * sigSt * (trap0 + trap1) / (6. * _PI_ * _PI_);
        }

        maxBndNu = 20. * tnuH;
        if (maxBndNu < (10 * Mnu)){maxBndNu=10*Mnu;}
        nNu = 0.;
        rhoNu = 0.;
        presNu = 0.;
        // Nu integration
        sigSt = (double)((maxBndNu - Mnu) / sigP);
        for (int ii=0; ii < (sigP-1); ii++) {
            ehold0 = Mnu + sigSt * (double)ii;
            ehold1 = Mnu + sigSt * (double)(ii+1);
            if (ii==0) {ehold0 += 1e-20;}
            // n
            trap0 =  ehold0 * sqrt(ehold0*ehold0 - Mnu*Mnu) /  ( exp((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (2 * _PI_ * _PI_);
            trap1 = ehold1 * sqrt(ehold1*ehold1 - Mnu*Mnu) /  ( exp((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (2 * _PI_ * _PI_);
            nNu += 0.5 * sigSt * (trap0 + trap1);
            // rho
            trap0 =  ehold0 * ehold0 * sqrt(ehold0*ehold0 - Mnu*Mnu) /  ( exp((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (2 * _PI_ * _PI_);
            trap1 =  ehold1 * ehold1 * sqrt(ehold1*ehold1 - Mnu*Mnu) /  ( exp((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (2 * _PI_ * _PI_);
            rhoNu += 0.5 * sigSt * (trap0 + trap1);
            // pressure
            trap0 =   pow(ehold0*ehold0 - Mnu*Mnu, 1.5) /  ( exp((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (6. * _PI_ * _PI_);
            trap1 =   pow(ehold1*ehold1 - Mnu*Mnu, 1.5) /  ( exp((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (6. * _PI_ * _PI_);
            presNu += 0.5 * sigSt * (trap0 + trap1);
        }

        neff = (degNu * rhoNu + rhoMaj) / h_cmb * (8./7.) * pow(11./4., 4./3.)  ;
        //printf("indx: %d \t z: %e \t Tcmb: %e \t RhoM/RhoN: %e \t Neff: %e \n", *lenIndx, zstore, TcmbStore, rhoMaj/ (degNu * rhoNu), neff );
        if ((rhoMaj / (degNu * rhoNu) < 1e-6)&&(TcmbStore < mMaj)){break;}
        if (GammaEff > 1){
            if (*lenIndx < 1){delT *= 5.;};}
        else {
        if (*lenIndx < 1){delT *= 4;};}


        pba->z_maj[*lenIndx] = zstore;
        pba->T_maj[*lenIndx] = TMstore;
        pba->T_nu[*lenIndx] = TNstore;
        pba->Mu_maj[*lenIndx] = muMstore;
        pba->Mu_nu[*lenIndx] = muNstore;
        //z_maj[*lenIndx] = zstore;
        //rho_maj[*lenIndx] = rhoMaj;
        //rho_nu[*lenIndx] = degNu * rhoNu;
        //press_maj[*lenIndx] = presMaj;
        //press_nu[*lenIndx] = degNu * presNu;

        *lenIndx+=1;
    }
    //printf("indx: %d \t z: %e \t Tcmb: %e \t Tnu: %e \t Tmaj: %e \t Munu: %e \t MuMaj: %e \n", indx, zstore,  TcmbStore, TNstore, TMstore, muNstore, muMstore);

    if (TcmbStore < (0.01 * mMaj) ) {
    //printf("Sucessful Finish... \n");
        stop_loop = true;}

  }

    *lenIndx -= 1;
//    z_maj[lenIndx] = zhold;
//    rho_maj[lenIndx] = rho_maj[lenIndx - 1] * 1e-3;
//    rho_nu[lenIndx] = rho_nu[lenIndx - 1] + rho_maj[lenIndx - 1];
//    press_maj[lenIndx] = press_maj[lenIndx - 1] * 1e-3;
//    press_nu[lenIndx] = press_nu[lenIndx - 1];

  return _SUCCESS_;
}

int RK_Eval(struct background *pba, double GammaPhi, double zhold, double tmajH, double tnuH, double muMh, double muNh, double mMaj, double tcur, double Mnu, double k[5]) {
    double epsil=1e-10;
    double H_preF = 8.*_PI_/3. * 6.707e-57; // eV^-2
    //double minMuChP = -1e-50;
    double nMaj=0, rhoMaj=0, presMaj=0, dnPdT=0, drPdT=0, dnPdmu=0, drPdmu=0, ColTrho=0, ColTn=0, sigSt;
    double ehold0, ehold1, trap0, trap1, dnphidt, drhophidt, dnNudt, drhoNudt;
    double nNu=0, rhoNu=0, presNu=0, dnNdT=0, drNdT=0, dnNdmu=0, drNdmu=0, maxBndNu, maxBndMaj;
    double h_cmb, h_mat, holdCT, degNu=6, degMaj=1, Hub;
    int sigP=100;
    bool mu_off= false;

//    if (fabs(muMh) < fabs(minMuChP)){muMh = minMuChP;}
//    if (!isfinite(muMh)||(muMh>=minMuChP)) {muMh = minMuChP;}
//    if (!isfinite(muNh)||(muNh>=minCHP)) {muNh = minCHP;}
//    if (!isfinite(tmajH)||(tmajH <= minMT)) {tmajH = minMT;}
//    checkRat = fabs(muMh / tmajH);
    //if (fabs(muMh / tmajH) > 40.) {muMh /= 2.;}


    h_cmb = pow(_PI_, 2.) / 15. * pow(tcur, 4.); // energy density in ev^4
    h_mat = (pba->Omega0_b + pba->Omega0_cdm) * pow((1. + zhold), 3) * 8.0835e-11 * pow(2.998e10 * 6.58e-16, 3) * pow(pba->h, 2); // eV^4

    maxBndMaj = 20. * tmajH;
    if (maxBndMaj < (3. * mMaj)){maxBndMaj=3. * mMaj;}
    nMaj = 0.;
    rhoMaj = 0.;
    presMaj = 0.;
    dnPdT=0.;
    drPdT=0.;
    dnPdmu=0.;
    drPdmu=0.;
    ColTrho=0.;
    ColTn=0.;
    // Maj integration
    sigSt = (double)((maxBndMaj - mMaj) / sigP);
    for (int ii=0; ii < (sigP-1); ii++) {
        ehold0 = mMaj + sigSt * (double)ii;
        ehold1 = mMaj + sigSt * (double)(ii+1);
        if (ii==0) {ehold0 += epsil;}
        // n
        trap0 = ehold0 * sqrt(ehold0*ehold0 - mMaj*mMaj) /  ( exp((ehold0 - muMh) / tmajH) - 1.) ;
        trap1 = ehold1 * sqrt(ehold1*ehold1 - mMaj*mMaj) /  ( exp((ehold1 - muMh) / tmajH) - 1.);
        nMaj += degMaj * 0.5 * sigSt * (trap0 + trap1)  / (2. * _PI_ * _PI_);
        // rho
        trap0 = ehold0 * ehold0 * sqrt(ehold0*ehold0 - mMaj*mMaj) /  ( exp((ehold0 - muMh) / tmajH) - 1.);
        trap1 = ehold1 * ehold1 * sqrt(ehold1*ehold1 - mMaj*mMaj) /  ( exp((ehold1 - muMh) / tmajH) - 1.) ;
        rhoMaj += degMaj * 0.5 * sigSt * (trap0 + trap1)  / (2. * _PI_ * _PI_);
        // Pres
        trap0 =   pow(ehold0*ehold0 - mMaj*mMaj, 1.5) /  ( exp((ehold0 - muMh) / tmajH) - 1.) ;
        trap1 =  pow(ehold1*ehold1 - mMaj*mMaj, 1.5) /  ( exp((ehold1 - muMh) / tmajH) - 1.) ;
        presMaj += degMaj * 0.5 * sigSt * (trap0 + trap1) / (6. * _PI_ * _PI_);
        // dndT
        trap0 = ehold0 * sqrt(ehold0*ehold0 - mMaj*mMaj) * ((ehold0 - muMh) / (4. * tmajH *tmajH)) /   pow(sinh((ehold0 - muMh) / 2. / tmajH), 2.) ;
        trap1 = ehold1 * sqrt(ehold1*ehold1 - mMaj*mMaj) * ((ehold1 - muMh) / (4. * tmajH *tmajH)) /   pow(sinh((ehold1 - muMh) / 2. / tmajH), 2.) ;
        dnPdT += degMaj * 0.5 * sigSt * (trap0 + trap1) / (2. * _PI_ * _PI_);
        // drhodT
        trap0 = ehold0 * ehold0 * sqrt(ehold0*ehold0 - mMaj*mMaj) * ((ehold0 - muMh) / (4. * tmajH *tmajH)) /   pow(sinh((ehold0 - muMh) / 2. / tmajH), 2.) ;
        trap1 = ehold1 * ehold1 * sqrt(ehold1*ehold1 - mMaj*mMaj) * ((ehold1 - muMh) / (4. * tmajH *tmajH)) /   pow(sinh((ehold1 - muMh) / 2. / tmajH), 2.) ;
        drPdT += degMaj * 0.5 * sigSt * (trap0 + trap1) / (2. * _PI_ * _PI_);
        // dndmu
        trap0 = ehold0 * sqrt(ehold0*ehold0 - mMaj*mMaj) * (1. / (4.*tmajH)) /   pow(sinh((ehold0 - muMh) / 2. / tmajH), 2.) ;
        trap1 = ehold1 * sqrt(ehold1*ehold1 - mMaj*mMaj) * (1. / (4. * tmajH )) /   pow(sinh((ehold1 - muMh) / 2. / tmajH), 2.) ;
        dnPdmu += degMaj * 0.5 * sigSt * (trap0 + trap1) / (2. * _PI_ * _PI_);
        // drdmu
        trap0 = ehold0 * ehold0 * sqrt(ehold0*ehold0 - mMaj*mMaj) * (1 / (4. * tmajH)) / pow(sinh((ehold0 - muMh) / 2. / tmajH), 2.) ;
        trap1 = ehold1 * ehold1 * sqrt(ehold1*ehold1 - mMaj*mMaj) * (1 / (4. * tmajH )) / pow(sinh((ehold1 - muMh) / 2. / tmajH), 2.) ;
        drPdmu += degMaj *0.5 * sigSt * (trap0 + trap1) / (2. * _PI_ * _PI_);
        // Collision term
        trap0 = -ehold0 * GammaPhi * mMaj / ehold0 * tnuH * (exp(ehold0 / tmajH + 2. *muNh/tnuH) - exp(ehold0 / tnuH + muMh/tmajH)) / ((exp(ehold0 / tnuH) - exp(2.*muNh/tnuH))*(exp(ehold0/tmajH) - exp(muMh/tmajH)));
        holdCT = ( exp((ehold0 - sqrt(ehold0*ehold0 - mMaj *mMaj)) / 2. / tnuH) + exp(muNh / tnuH)) *  (exp((ehold0 + sqrt(ehold0*ehold0 - mMaj *mMaj) + 2. * muNh) / 2. / tnuH) + exp(ehold0 / tnuH));
        holdCT /= ((exp((ehold0 + sqrt(ehold0*ehold0 - mMaj *mMaj)) / 2. / tnuH) + exp(muNh / tnuH) ) * (exp((ehold0 + sqrt(ehold0*ehold0 - mMaj *mMaj) + 2. * muNh) / 2. / tnuH) + exp(ehold0 / tnuH)));
        trap0 *= log(holdCT);

        trap1 = -ehold1 * GammaPhi * mMaj / ehold1 * tnuH * (exp(ehold1 / tmajH + 2. *muNh/tnuH) - exp(ehold1 / tnuH + muMh/tmajH)) / ((exp(ehold1 / tnuH) - exp(2.*muNh/tnuH))*(exp(ehold1/tmajH) - exp(muMh/tmajH)));
        holdCT = ( exp((ehold1 - sqrt(ehold1*ehold1 - mMaj *mMaj)) / 2. / tnuH) + exp(muNh / tnuH)) *  (exp((ehold1 + sqrt(ehold1*ehold1 - mMaj *mMaj) + 2. * muNh) / 2. / tnuH) + exp(ehold1 / tnuH));
        holdCT /= ((exp((ehold1 + sqrt(ehold1*ehold1 - mMaj *mMaj)) / 2. / tnuH) + exp(muNh / tnuH) ) * (exp((ehold1 + sqrt(ehold1*ehold1 - mMaj *mMaj) + 2. * muNh) / 2. / tnuH) + exp(ehold1 / tnuH)));
        trap1 *= log(holdCT);
//        trap0 = - GammaPhi * sqrt(ehold0*ehold0 - mMaj*mMaj) * mMaj * (exp((muMh - ehold0)/tmajH) - exp((2*muNh - ehold0)/tnuH));
//        trap1 = - GammaPhi * sqrt(ehold1*ehold1 - mMaj*mMaj) * mMaj * (exp((muMh - ehold1)/tmajH) - exp((2*muNh - ehold1)/tnuH));

        ColTn += 0.5 * sigSt * (trap0 + trap1) / (2. * _PI_ * _PI_);
        ColTrho +=  0.5 * sigSt * (trap0*ehold0 + trap1*ehold1) / (2. * _PI_ * _PI_);
    }
    dnphidt = ColTn;
    drhophidt = ColTrho;
    if (isnan(dnphidt)){dnphidt = 0.;}
    if (isnan(drhophidt)){drhophidt = 0.;}
    if (!isfinite(rhoMaj)){rhoMaj = 0.;}
    if (!isfinite(nMaj)){nMaj = 0.;}

    drhoNudt = -drhophidt / 6.;
    dnNudt = -2.*dnphidt / 6.;

    maxBndNu = 20. * tnuH;
    if (maxBndNu < (5. * Mnu)){maxBndNu=5.*Mnu;}
    // Nu integration
    sigSt = (double)((maxBndNu - Mnu) / sigP);
    for (int ii=0; ii < (sigP-1); ii++) {
        ehold0 = Mnu + sigSt * (double)ii;
        ehold1 = Mnu + sigSt * (double)(ii+1);
        if (ii==0) {ehold0 += epsil;}
        // n
        trap0 =   ehold0 * sqrt(ehold0*ehold0 - Mnu*Mnu) /  ( exp((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (2. * _PI_ * _PI_);
        trap1 = ehold1 * sqrt(ehold1*ehold1 - Mnu*Mnu) /  ( exp((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (2. * _PI_ * _PI_);
        nNu += 0.5 * sigSt * (trap0 + trap1);
        // rho
        trap0 =   ehold0 * ehold0 * sqrt(ehold0*ehold0 - Mnu*Mnu) /  ( exp((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (2. * _PI_ * _PI_);
        trap1 =  ehold1 * ehold1 * sqrt(ehold1*ehold1 - Mnu*Mnu) /  ( exp((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (2. * _PI_ * _PI_);
        rhoNu += 0.5 * sigSt * (trap0 + trap1);
        // pressure
        trap0 =  pow(ehold0*ehold0 - Mnu*Mnu, 1.5) /  ( exp((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (6. * _PI_ * _PI_);
        trap1 =  pow(ehold1*ehold1 - Mnu*Mnu, 1.5) /  ( exp((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / tnuH) + 1.) / (6. * _PI_ * _PI_);
        presNu += 0.5 * sigSt * (trap0 + trap1);
        // dndT
        trap0 =  ehold0 * sqrt(ehold0*ehold0 - Mnu*Mnu) * ((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / (4 * tnuH *tnuH)) /   pow(cosh((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / 2. / tnuH), 2.) ;
        trap1 =  ehold1 * sqrt(ehold1*ehold1 - Mnu*Mnu) * ((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / (4 * tnuH *tnuH)) /   pow(cosh((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / 2. / tnuH), 2.) ;
        dnNdT += 0.5 * sigSt * (trap0 + trap1) / (2. * _PI_ * _PI_);
        // drhodT
        trap0 = ehold0 * ehold0 * sqrt(ehold0*ehold0 - Mnu*Mnu) * ((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / (4. * tnuH *tnuH)) /   pow(cosh((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / 2. / tnuH), 2.) ;
        trap1 =  ehold1 * ehold1 * sqrt(ehold1*ehold1 - Mnu*Mnu) * ((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / (4. * tnuH *tnuH)) /   pow(cosh((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / 2. / tnuH), 2.) ;
        drNdT += 0.5 * sigSt * (trap0 + trap1) / (2. * _PI_ * _PI_);
        // dndmu
        trap0 = ehold0 * sqrt(ehold0*ehold0 - Mnu*Mnu)  /   (2. * tnuH * cosh((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / tnuH) + 2. *tnuH) ;
        trap1 = ehold1 * sqrt(ehold1*ehold1 - Mnu*Mnu) /   (2. * tnuH * cosh((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / tnuH) + 2. *tnuH) ;
        dnNdmu += 0.5 * sigSt * (trap0 + trap1) / (2. * _PI_ * _PI_);
        // drdmu
        trap0 =  ehold0 * ehold0 * sqrt(ehold0*ehold0 - Mnu*Mnu) /  (2. * tnuH * cosh((sqrt(ehold0*ehold0 - Mnu*Mnu) - muNh) / tnuH) + 2. *tnuH) ;
        trap1 =   ehold1 * ehold1 * sqrt(ehold1*ehold1 - Mnu*Mnu) / (2. * tnuH * cosh((sqrt(ehold1*ehold1 - Mnu*Mnu) - muNh) / tnuH) + 2. *tnuH) ;
        drNdmu += 0.5 * sigSt * (trap0 + trap1) / (2. * _PI_ * _PI_);
    }
    Hub = sqrt( H_preF * (degNu * rhoNu + rhoMaj + h_mat + h_cmb)) ; // eV

    k[2] = -(3*((presMaj + rhoMaj) * dnPdT - nMaj * drPdT) - (dnPdT * drhophidt - drPdT * dnphidt) / Hub) / (dnPdmu * drPdT - dnPdT * drPdmu) / (1. + zhold) / muMh;
    k[0] = (3*((presMaj + rhoMaj) * dnPdmu - nMaj * drPdmu) - (dnPdmu * drhophidt - drPdmu * dnphidt) / Hub) / (dnPdmu * drPdT - dnPdT * drPdmu) / (1. + zhold) / tmajH;
    k[1] = (3*((presNu + rhoNu) * dnNdmu - nNu * drNdmu) - (dnNdmu * drhoNudt - drNdmu * dnNudt) / Hub) / (dnNdmu * drNdT - dnNdT * drNdmu) / (1. + zhold) / tnuH;
    k[3] = -(3*((presNu + rhoNu) * dnNdT - nNu * drNdT) - (dnNdT * drhoNudt - drNdT * dnNudt) / Hub) / (dnNdmu * drNdT - dnNdT * drNdmu) / (1. + zhold) / muNh;
    k[4] = tcur / (1. + zhold);
    //printf("Check: %e \t %e \t %e \t %e \t %e \n", tmajH, muMh, rhoMaj, nMaj,presMaj);

    return _SUCCESS_; }

double bessk( int n, double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n >= 0*/
/* Note that for x == 0 the functions bessy and bessk are not */
/* defined and a blank is returned.                           */
/*------------------------------------------------------------*/
{
   int j;
   double bk,bkm,bkp,tox;


   if (n < 0 || x == 0.0)
   {
      return( 0. );
   }
   if (n == 0)
      return( bessk0(x) );
   if (n == 1)
      return( bessk1(x) );

   tox=2.0/x;
   bkm=bessk0(x);
   bk=bessk1(x);
   for (j=1;j<n;j++) {
      bkp=bkm+j*tox*bk;
      bkm=bk;
      bk=bkp;
   }
   return bk;
}

double bessk0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=0.  */
/*------------------------------------------------------------*/
{
   double y,ans;

   if (x <= 2.0) {
      y=x*x/4.0;
      ans=(-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
         +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
         +y*(0.10750e-3+y*0.74e-5))))));
   } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
         +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
         +y*(-0.251540e-2+y*0.53208e-3))))));
   }
   return ans;
}




double bessk1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function Kn(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double y,ans;

   if (x <= 2.0) {
      y=x*x/4.0;
      ans=(log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
         +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
         +y*(-0.110404e-2+y*(-0.4686e-4)))))));
   } else {
      y=2.0/x;
      ans=(exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
         +y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
         +y*(0.325614e-2+y*(-0.68245e-3)))))));
   }
   return ans;
}


double bessi0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=0.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10


double bessi1( double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) and n=1.  */
/*------------------------------------------------------------*/
{
   double ax,ans;
   double y;


   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}

double bessj0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          0 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         -y*0.934935152e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
   }
   return ans;
}



double bessj1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          1 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
         +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
      ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
         +y*(99447.43394+y*(376.9991397+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-2.356194491;
      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
      ans2=0.04687499995+y*(-0.2002690873e-3
         +y*(0.8449199096e-5+y*(-0.88228987e-6
         +y*0.105787412e-6)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
      if (x < 0.0) ans = -ans;
   }
   return ans;
}



/*
#>            bessj.dc2

Function:     bessj

Purpose:      Evaluate Bessel function of first kind of integer order.

Category:     MATH

File:         bessel.c

Author:       M.G.R. Vogelaar

Use:          #include "bessel.h"
              double   result;
              result = bessj( int n,
                              double x )


              bessj    Return the Bessel function of integer order
                       for input value x.
              n        Integer order of Bessel function.
              x        Double at which the function is evaluated.


Description:  bessj evaluates at x the Bessel function of the first kind
              and of integer order n.
              This routine is NOT callable in FORTRAN.

Updates:      Jun 29, 1998: VOG, Document created.
#<
*/


double bessj( int n, double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          n at input x                                      */
/* The function can also be called for n = 0 and n = 1.       */
/*------------------------------------------------------------*/
{
   int    j, jsum, m;
   double ax, bj, bjm, bjp, sum, tox, ans;


   if (n < 0)
   {
      double   dblank;
      return( 0. );
   }
   ax=fabs(x);
   if (n == 0)
      return( bessj0(ax) );
   if (n == 1)
      return( bessj1(ax) );


   if (ax == 0.0)
      return 0.0;
   else if (ax > (double) n) {
      tox=2.0/ax;
      bjm=bessj0(ax);
      bj=bessj1(ax);
      for (j=1;j<n;j++) {
         bjp=j*tox*bj-bjm;
         bjm=bj;
         bj=bjp;
      }
      ans=bj;
   } else {
      tox=2.0/ax;
      m=2*((n+(int) sqrt(ACC*n))/2);
      jsum=0;
      bjp=ans=sum=0.0;
      bj=1.0;
      for (j=m;j>0;j--) {
         bjm=j*tox*bj-bjp;
         bjp=bj;
         bj=bjm;
         if (fabs(bj) > BIGNO) {
            bj *= BIGNI;
            bjp *= BIGNI;
            ans *= BIGNI;
            sum *= BIGNI;
         }
         if (jsum) sum += bj;
         jsum=!jsum;
         if (j == n) ans=bjp;
      }
      sum=2.0*sum-bj;
      ans /= sum;
   }
   return  x < 0.0 && n%2 == 1 ? -ans : ans;
}

double bessi( int n, double x)
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate modified Bessel function In(x) for n >= 0*/
/*------------------------------------------------------------*/
{
   int j;
   double bi,bim,bip,tox,ans;


   if (n < 0)
   {
      double   dblank;

      return( 0. );
   }
   if (n == 0)
      return( bessi0(x) );
   if (n == 1)
      return( bessi1(x) );


   if (x == 0.0)
      return 0.0;
   else {
      tox=2.0/fabs(x);
      bip=ans=0.0;
      bi=1.0;
      for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
         bim=bip+j*tox*bi;
         bip=bi;
         bi=bim;
         if (fabs(bi) > BIGNO) {
            ans *= BIGNI;
            bi *= BIGNI;
            bip *= BIGNI;
         }
         if (j == n) ans=bip;
      }
      ans *= bessi0(x)/bi;
      return  x < 0.0 && n%2 == 1 ? -ans : ans;
   }
}
