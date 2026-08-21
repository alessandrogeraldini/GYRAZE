/*
 * densfinorb_par.c — OpenMP-parallel version of densfinorb.
 *
 * Identical to densfinorb except the closed-orbit j (xbar) loop is
 * parallelised via OpenMP.  The six per-j accumulators that drove the
 * sequential dependency (intdvx, intdvx_corr_delta, intdvx_corr_chiM,
 * intdvx_in, intdvx_ref) are stored per-j, filled in parallel, then
 * reduced sequentially with the trapezoidal rule.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include "mps.h"

#define APPROXMUFORSMALLPHI 0
#define TESTELL 0
#ifndef numb
#define numb 0.00000001
#endif

/* -----------------------------------------------------------------------
 * Static helpers (copies of the same functions in denscalc.c).
 * ----------------------------------------------------------------------- */
static double uperp_from_mu2(int j, double mu_target,
                              double **mu, double **Uperp, int lowerlimit, int upperlimit, int maxk,
                              int current_j, int current_k, int *k_out, double phi_imin)
{
    int k;
    (void)current_j; (void)current_k; (void)phi_imin;
    if (upperlimit < 0) return -1.0;
    double mu_min = 1e20, Uperp_min = 1e20, mu_max = -1.0;
    for (k = lowerlimit; k <= maxk; k++) {
        if (mu[j][k] < mu_min) { mu_min = mu[j][k]; Uperp_min = Uperp[j][k]; }
        if (mu[j][k] > mu_max)   mu_max = mu[j][k];
        double lo = mu[j][k+1], hi = mu[j][k];
        if (lo > hi) { double tmp = lo; lo = hi; hi = tmp; }
        if (lo <= mu_target && mu_target <= hi) {
            double dmu = mu[j][k+1] - mu[j][k];
            if (k_out) *k_out = k;
            if (dmu == 0.0) return Uperp[j][k];
            double t = (mu_target - mu[j][k]) / dmu;
            return Uperp[j][k] + t * (Uperp[j][k+1] - Uperp[j][k]);
        }
    }
    (void)mu_min; (void)Uperp_min; (void)mu_max;
    return -1234;
}

static double uperp_from_mu3(int j, double mu_target, double *xbarr, double *xx, double *phi, int size_xx)
{
    int i, i_closest = 0;
    double best = fabs(xx[0] - xbarr[j]);
    for (i = 1; i < size_xx; i++) {
        double d = fabs(xx[i] - xbarr[j]);
        if (d < best) { best = d; i_closest = i; }
    }
    double phi_pp;
    if (i_closest == 0) {
        double dx = xx[1] - xx[0];
        phi_pp = (phi[2] - 2*phi[1] + phi[0]) / (dx*dx);
    } else if (i_closest == size_xx - 1) {
        double dx = xx[size_xx-1] - xx[size_xx-2];
        phi_pp = (phi[size_xx-1] - 2*phi[size_xx-2] + phi[size_xx-3]) / (dx*dx);
    } else {
        double dxl = xx[i_closest] - xx[i_closest-1];
        double dxr = xx[i_closest+1] - xx[i_closest];
        phi_pp = 2.0*(phi[i_closest+1]/dxr - phi[i_closest]*(1.0/dxl + 1.0/dxr) + phi[i_closest-1]/dxl) / (dxl + dxr);
    }
    return mu_target * sqrtf(1.0 + phi_pp) + phi[i_closest];
}

/* -----------------------------------------------------------------------
 * densfinorb_par
 * ----------------------------------------------------------------------- */
void densfinorb_par(double Ti, double lenfactor, double alpha,
                    int size_phigrid, int *size_ngrid,
                    double *n_grid, double *n_grid_corr_delta,
                    double *n_grid_corr_chiM,
                    double *x_grid, double *phi_grid, double charge,
                    double **FF, double *mumu, double *UU,
                    int sizemumu, int sizeUU,
                    double grid_parameter, double *flux, double *Qflux,
                    int zoomfactor, double margin, double phi_DSbump,
                    double *vy_op, double *mu_op, double *chiMax_op,
                    double *dmudvy_op, int *size_op,
                    FILE *fmu, FILE *fjmc_out)
{
    /* ================================================================
     * PHASE 1 — ARRAY FILLING  (identical to densfinorb)
     * ================================================================ */
    clock_t begin = clock();
    double wt_begin = omp_get_wtime();
    double limit_rho = 8.0;
    double n_inf = 0.0;
    double Ucrit = 0.0, Ucritold = 0.0, Ucritp = 0.0, Ucritpold = 0.0, frac_reflected = 0.0, frac = 0.0, vzcrit;
    double deltax, deltax_inf, deltaE = 0.1, phibar;
    double *chiinf, *muinf, **vxinf;
    double *xx, *phi, *phip, *phipp, **chi;
    int s, w, ind, muapproxtype = 1;
    int stop = 0, sizexbar, maxj, j_inf;
    int reflected = 0;
    int icrit;
    double *openorbit, openorbitantycal, **mu, *muopen, *chiMopen, *xbaropen, *openorbitopen, *Ucritf, *xbar, xbarcrit, chiMcrit, *xifunction;
    int i = 0, ic = 0, j = 0, k = 0, l = 0;
    int phi_monotone = 1, phi_imax = -1, phi_imin = -1;
    int *jmclosed, *jmopen, sizeU, size_finegrid, size_xlim;
    int *crossed_max, *crossed_min, *kdrop;
    int *lowerlimit, *upperlimit, **upper, *imax, *imin;
    double **Uperp, ***vx, *chiMax, *chimpp, *chimin, oorbintgrd, oorbintgrdantycal;
    double vz, U, dvz = 0.1, dvzopen = 0.1, dvx, dxbar, intdU = 0.0, intdUopen = 0.0, intdU_corr_delta = 0.0, intdU_corr_chiM = 0.0;
    double intdUold = 0.0, intdU_corr_delta_old = 0.0, intdU_corr_chiM_old = 0.0, intdvx = 0.0, intdvxold = 0.0, intdvx_corr_delta = 0.0, intdvx_corr_delta_old = 0.0, intdvx_corr_chiM = 0.0, intdvx_corr_chiM_old = 0.0, intdxbar = 0.0, intdxbar_corr_delta = 0.0, intdxbar_corr_chiM = 0.0, intdxbaropen = 0.0, F, Fold = 0.0, Fold_ref = 0.0, Ucap;
    double intdUopenflow = 0.0, intdUopenflowold = 0.0, intdxbaropenflow = 0.0, oorbintgrdflow = 0.0, oorbintgrdflowold = 0.0, oorbintgrdener = 0.0, oorbintgrdenerold = 0.0;
    double intdU_in = 0.0, intdU_in_old = 0.0, intdU_ref = 0.0, intdU_ref_old = 0.0;
    double intdvx_in = 0.0, intdvx_in_old = 0.0, intdvx_ref = 0.0, intdvx_ref_old = 0.0;
    double intdxbar_in = 0.0, intdxbar_ref = 0.0;
    FILE *fout_inref = NULL;
    double oorbintgrdold = 0.0, Fopen = 0.0, intdUopenold = 0.0, intdUopenener = 0.0, intdUopenenerold = 0.0, intdxbaropenener = 0.0, Qflux0 = 0.0;
    double vx0open;
    double intdUantycal = 0.0, intdvxantycal = 0.0, vxnew = 0.0, vxold = 0.0, Uperpnew = 0.0, U_lb = 0.0, *xtop, intdUopenantycal = 0.0;
    double openorbitnew, chinew, munew = 0.0, muold = 0.0, Uperpold = 0.0, dmu_dUperp = 0.0, dmu_dUperp_old = 0.0;
    double xi, *gg, *ff;
    double flux0, du, fluxinf1old, fluxinfintgrdold, fluxinfintgrd, Qfluxinf1old, Qfluxinfintgrdold, Qfluxinfintgrd, u, Chodura2, Chodura2old, Chodura1old, Chodura1, Chodura;
    double fluxinf, Qfluxinf, fluxinf1, Qfluxinf1, densinf1, densinf, densinf1old;
    double musmall, muell, Omegaell, minphiformucalc = 0.3;
    FILE *fp = NULL, *filellip = NULL;

    gg = malloc(size_phigrid * sizeof(double));
    ff = malloc(size_phigrid * sizeof(double));
    for (i = 0; i < size_phigrid; i++) {
        gg[i] = sqrt(x_grid[i]);
        ff[i] = pow(sqrt(grid_parameter) + gg[i], 2.0) - grid_parameter;
    }
    deltax = ff[1];
    size_finegrid = zoomfactor * size_phigrid - zoomfactor;
    printf("size_finegrid = %d\t it is %d x %d - %d\n", size_finegrid, zoomfactor, size_phigrid, zoomfactor);
    xx  = (double*)calloc(size_finegrid, sizeof(double));
    phi = (double*)calloc(size_finegrid, sizeof(double));
    filellip = fopen("checkellip.txt", "w");
    if (filellip == NULL) printf("Cannot open checkellip.txt\n");
    if (zoomfactor != 1) {
        i = 0;
        gsl_interp_accel *acc = gsl_interp_accel_alloc();
        gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, size_phigrid);
        gsl_spline_init(spline, ff, phi_grid, size_phigrid);
        fp = fopen("OUTPUT/phispline.txt", "w");
        if (fp == NULL) printf("Error: phispline not created\n");
        for (i = 0; i < size_finegrid; i++) {
            xi = i * deltax / zoomfactor;
            if (i == 0) xi += 0.00001;
            if (i == size_finegrid - 1) xi -= 0.00001;
            xx[i]  = pow(pow(grid_parameter + xi, 0.5) - sqrt(grid_parameter), 2.0);
            phi[i] = gsl_spline_eval(spline, xi, acc);
            if (fp != NULL) fprintf(fp, "%f %f\n", xx[i], phi[i]);
            xx[i]  *= lenfactor;
            phi[i] *= (charge / Ti);
        }
        if (fp != NULL) fclose(fp);
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    } else {
        fp = fopen("OUTPUT/phispline.txt", "w");
        if (fp == NULL) printf("Error: phispline not created\n");
        for (i = 0; i < size_finegrid; i++) {
            xi    = i * deltax / zoomfactor;
            xx[i] = x_grid[i] * lenfactor;
            phi[i] = phi_grid[i];
            if (fp != NULL) fprintf(fp, "%f %f\n", xx[i], phi[i]);
            phi[i] *= (charge / Ti);
        }
        if (fp != NULL) fclose(fp);
    }
    printf("lenfactor = %f\n", lenfactor);
    printf("Ti= %f\tcharge = %f\n", Ti, charge);
    printf("phi[0] = %f\tphi_grid[0] = %f\n", phi[0], phi_grid[0]);

    if (charge < 0) {
        double dphi0 = phi[1] - phi[0];
        double dphi_prev = dphi0;
        for (i = 1; i < size_finegrid - 1; i++) {
            double dphi = phi[i+1] - phi[i];
            if (dphi * dphi_prev < 0.0) {
                phi_monotone = 0;
                if (dphi < 0.0) phi_imin = i;
                else             phi_imax = i;
            }
            dphi_prev = dphi;
        }
        if (phi_monotone != 1)
            printf("WARNING in densfinorb_par: phi is not monotone "
                   "(local max at i=%d, x=%.6f; local min at i=%d, x=%.6f)\n",
                   phi_imax, xx[phi_imax], phi_imin, xx[phi_imin]);
        else
            printf("phi monotonicity check passed\n");
    }

    Ucap = 12.0 + 10.0 / Ti;
    if (bilin_interp(0.0, Ucap, FF, mumu, UU, sizemumu, sizeUU, -1, -1) > 1e-6)
        printf("ERROR in densfinorb_par: increase Ucap please\n");
    else if (bilin_interp(Ucap, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1) > 1e-6)
        printf("ERROR in densfinorb_par: increase Ucap please\n");
    else
        printf("Ucap has been checked to be large enough\n");

    Qfluxinfintgrd = Qfluxinf1 = Qfluxinf = fluxinfintgrd = fluxinf1 = fluxinf = densinf = densinf1 = Chodura2 = Chodura1 = Chodura = 0.0;
    du = 0.01;
    printf("sizemumu = %d\tsizeUU=%d\n", sizemumu, sizeUU);
    for (i = 0; i < sizemumu; i++) {
        fluxinf1old = fluxinf1;
        Qfluxinf1old = Qfluxinf1;
        Chodura1old = Chodura1;
        Chodura1 = 0.0;
        fluxinf1 = 0.0;
        Qfluxinf1 = 0.0;
        densinf1old = densinf1;
        densinf1 = 0.0;
        Fold = F = 0.0;
        fluxinfintgrdold = fluxinfintgrd = 0.0;
        for (j = 0; j < sizeUU; j++) {
            u = sqrt(2.0 * UU[j]);
            fluxinfintgrdold  = fluxinfintgrd;
            Qfluxinfintgrdold = Qfluxinfintgrd;
            Chodura2old = Chodura2;
            Fold = F;
            F = FF[i][j];
            fluxinfintgrd  = F * u;
            Qfluxinfintgrd = F * (mumu[i] + 0.5*u*u) * u;
            if (j > 0) Chodura2 = F / (u*u);
            else       Chodura2 = 0.0;
            if (j != 0) {
                du = u - sqrt(2.0 * UU[j-1]);
                fluxinf1  += 0.5*(fluxinfintgrd  + fluxinfintgrdold)  * du;
                Qfluxinf1 += 0.5*(Qfluxinfintgrd + Qfluxinfintgrdold) * du;
                Chodura1  += 0.5*(Chodura2        + Chodura2old)       * du;
                densinf1  += 0.5*(F               + Fold)              * du;
            }
        }
        if (i != 0) {
            Qfluxinf += 0.5*(mumu[i] - mumu[i-1]) * (Qfluxinf1 + Qfluxinf1old);
            fluxinf  += 0.5*(mumu[i] - mumu[i-1]) * (fluxinf1  + fluxinf1old);
            Chodura  += 0.5*(mumu[i] - mumu[i-1]) * (Chodura1  + Chodura1old);
            densinf  += 0.5*(mumu[i] - mumu[i-1]) * (densinf1  + densinf1old);
        }
    }
    densinf   *= (4.0 * M_PI);
    fluxinf   *= (4.0 * M_PI);
    Qfluxinf  *= (4.0 * M_PI);
    Chodura   *= (4.0 * M_PI);
    *flux  = fluxinf  / densinf;
    *Qflux = Qfluxinf / densinf;
    Chodura /= densinf;
    printf("densinf = %f\tfluxinf = %f\tQfluxinf = %f\tChodura = %f\n", densinf, *flux, *Qflux, Chodura);
    printf("sizemumu = %d\tsizeUU = %d\n", sizemumu, sizeUU);

    phip       = malloc(size_finegrid * sizeof(double));
    phipp      = malloc(size_finegrid * sizeof(double));
    jmclosed   = malloc(size_finegrid * sizeof(int));
    jmopen     = malloc(size_finegrid * sizeof(int));
    xifunction = malloc(size_finegrid * sizeof(double));
    xbar = (double*)calloc(size_finegrid, sizeof(double));
    for (i = 0; i < size_finegrid; i++) {
        jmopen[i] = jmclosed[i] = 0;
        if (i == 0) {
            phip[0] = (phi[1] - phi[0]) / (xx[1] - xx[0]);
        } else if (i == size_finegrid - 1) {
            phip[i] = phip[i-1];
        } else {
            phip[i] = ((xx[i] - xx[i-1]) / (xx[i+1] - xx[i-1])) * (phi[i+1] - phi[i]) / (xx[i+1] - xx[i])
                    + ((xx[i+1] - xx[i]) / (xx[i+1] - xx[i-1])) * (phi[i]   - phi[i-1]) / (xx[i]   - xx[i-1]);
        }
        xifunction[i] = xx[i] + phip[i];
        if (DEBUG == 1) printf("xifunction[%d] = %f\n", i, xifunction[i]);
        if (i == 1) {
            if (xifunction[i] > xifunction[i-1]) icrit = i - 1;
        } else if (i > 1) {
            if ((xifunction[i] > xifunction[i-1]) && (xifunction[i-1] < xifunction[i-2]))
                icrit = i - 1;
            else if ((xifunction[i] < xifunction[i-1]) && (xifunction[i-1] > xifunction[i-2])) {
                printf("ERROR in densfinorb_par: too much noise in second derivative\n");
                printf("xifunction = %f\n", xifunction[i-1]);
                printf("i = %d\n", i);
            }
        }
    }

    for (i = 0; i < size_finegrid; i++) {
        if (i == 0)
            phipp[0] = (phip[1] - phip[0]) / (xx[1] - xx[0]);
        else if (i == size_finegrid - 1)
            phipp[i] = phipp[i-1];
        else
            phipp[i] = ((xx[i]   - xx[i-1]) / (xx[i+1] - xx[i-1])) * (phip[i+1] - phip[i])   / (xx[i+1] - xx[i])
                     + ((xx[i+1] - xx[i])   / (xx[i+1] - xx[i-1])) * (phip[i]   - phip[i-1]) / (xx[i]   - xx[i-1]);
    }

    j = 0;
    for (k = icrit + 1; k < size_finegrid; k++) {
        xbar[j] = xifunction[k];
        j++;
    }
    xbarcrit = xifunction[icrit];
    chiMcrit  = 0.5 * phip[icrit] * phip[icrit] + phi[icrit];
    sizexbar  = j;
    if (DEBUG == 1) for (j = 0; j < sizexbar; j++) printf("xbar[%d] = %f\n", j, xbar[j]);

    chi       = (double**)calloc(sizexbar, sizeof(double*));
    Uperp     = (double**)calloc(sizexbar, sizeof(double*));
    mu        = (double**)calloc(sizexbar, sizeof(double*));
    vx        = (double***)calloc(sizexbar, sizeof(double**));
    upper     = (int**)calloc(sizexbar, sizeof(int*));
    chimin    = (double*)calloc(sizexbar, sizeof(double));
    chiMax    = (double*)calloc(sizexbar, sizeof(double));
    chimpp    = (double*)calloc(sizexbar, sizeof(double));
    crossed_min = (int*)calloc(sizexbar, sizeof(int));
    crossed_max = (int*)calloc(sizexbar, sizeof(int));
    kdrop     = (int*)calloc(sizexbar, sizeof(int));
    openorbit = (double*)calloc(sizexbar, sizeof(double));
    upperlimit = (int*)calloc(sizexbar, sizeof(int));
    lowerlimit = (int*)calloc(sizexbar, sizeof(int));
    xtop      = (double*)calloc(sizexbar, sizeof(double));
    imax      = (int*)calloc(sizexbar, sizeof(int));
    imin      = (int*)calloc(sizexbar, sizeof(int));

    for (j = 0; j < sizexbar; j++) {
        imax[j] = imin[j] = -1;
        openorbit[j] = 0.0;
        crossed_min[j] = crossed_max[j] = 0;
        xtop[j] = chimin[j] = chiMax[j] = 0.0;
        upperlimit[j] = -1;
        lowerlimit[j] = 0;
    }
    for (j = 0; j < sizexbar; j++) {
        chi[j]   = (double*)calloc(size_finegrid, sizeof(double));
        Uperp[j] = (double*)calloc(size_finegrid, sizeof(double));
        mu[j]    = (double*)calloc(size_finegrid, sizeof(double));
        upper[j] = (int*)calloc(size_finegrid, sizeof(int));
        vx[j]    = (double**)calloc(size_finegrid, sizeof(double*));
        for (i = 0; i < size_finegrid; i++)
            vx[j][i] = (double*)calloc(size_finegrid, sizeof(double));
    }
    for (j = 0; j < sizexbar; j++)
        chi[j][0] = 0.5 * pow(xx[0] - xbar[j], 2.0) + phi[0];

    for (i = 1; i < size_finegrid; i++) {
        for (j = 0; j < sizexbar; j++) {
            chi[j][i] = 0.5 * pow(xx[i] - xbar[j], 2.0) + phi[i];
            if ((i == 1) && (chi[j][i] < chi[j][i-1])) {
                crossed_max[j] += 1;
                imax[j] = 0;
                if (phi_DSbump > 1e-6) {
                    chiMax[j] = chi[j][0] + phi_DSbump;
                    kdrop[j]  = ceil(phi_DSbump / deltaE);
                    deltaE    = phi_DSbump / kdrop[j];
                } else {
                    chiMax[j] = chi[j][0];
                    kdrop[j]  = 0;
                    deltaE    = 0.0;
                }
                for (k = 0; k <= kdrop[j]; k++) {
                    Uperp[j][k] = chiMax[j] - k * deltaE;
                    mu[j][k]    = 0.0;
                    vx[j][0][k] = sqrt(2.0 * Uperp[j][k] - chi[j][0]);
                }
                upper[j][0] = kdrop[j];
            } else if ((i > 1) && (chi[j][i] < chi[j][i-1]) && (chi[j][i-1] > chi[j][i-2])) {
                crossed_max[j] += 1;
                if (crossed_max[j] > 1) {
                    printf("***WARNING*** There is more than one maximum!\n");
                    printf("j is %d and maxima at %d and %d for first and second\n", j, imax[j], i-1);
                    for (k = imax[j]-1; k < i+1; k++) printf("chi[%d][%d] = %f\n", j, k, chi[j][k]);
                    crossed_max[j] = 1;
                    printf("exit code now\n");
                    exit(-1);
                }
                imax[j]   = i - 1;
                chiMax[j] = chi[j][i-1];
            } else if ((i > 1) && (chi[j][i] > chi[j][i-1]) && (chi[j][i-1] < chi[j][i-2])) {
                crossed_min[j] += 1;
                if (crossed_min[j] > 1) {
                    printf("***WARNING*** There is more than one minimum!\n");
                    printf("j is %d and minima at %d and %d for first and second\n", j, imin[j], i-1);
                    for (k = imin[j]-1; k < i+1; k++) printf("chi[%d][%d] = %f\n", j, k, chi[j][k]);
                    crossed_min[j] = 1;
                }
                imin[j]       = i - 1;
                upperlimit[j] = imin[j] - imax[j] + kdrop[j];
                chimpp[j]     = ((chi[j][i] - chi[j][i-1]) / (xx[i] - xx[i-1]) - (chi[j][i-1] - chi[j][i-2]) / (xx[i-1] - xx[i-2])) * 2.0 / (xx[i] - xx[i-2]);
                chimin[j]     = chi[j][i-1];
                crossed_max[j] = 1;
            }
            if (((crossed_max[j] == 1 && crossed_min[j] == 0) || (i-1 == imin[j]))) {
                Uperp[j][i-1-imax[j]+kdrop[j]] = chi[j][i-1];
                if (Uperp[j][i-1-imax[j]] > Ucap && lowerlimit[j] != 0)
                    lowerlimit[j] = i-1-imax[j]+kdrop[j];
                mu[j][i-1-imax[j]+kdrop[j]] = 0.0;
                upper[j][i-1] = i-1-imax[j]+kdrop[j];
                for (k = 0; k <= upper[j][i-1]; k++) {
                    if ((upper[j][i-1] == 0) || (i-1 == 0)) {
                        vx[j][i-1][k] = sqrt(2.0*Uperp[j][k] - chi[j][i-1]);
                        mu[j][k] += 0.0;
                    } else if ((k == imin[j] - imax[j] + kdrop[j]) && (crossed_min[j] == 1)) {
                        vx[j][i-1][k] = 0.0;
                        mu[j][k] = 0.0;
                    } else if ((k == imin[j] - imax[j] - 1 + kdrop[j]) && (crossed_min[j] == 1)) {
                        vx[j][i-1][k] = sqrt(2.0*(Uperp[j][k] - chi[j][i-1]));
                        mu[j][k] = 0.5 * pow(xx[imin[j]-1] - xx[imin[j]], 2.0) * pow(chimpp[j], 0.5);
                    } else if (k == upper[j][i-1] - 1) {
                        vx[j][i-1][k] = sqrt(2.0*(Uperp[j][k] - chi[j][i-1]));
                        mu[j][k] += (sqrt(2.0)/M_PI) * sqrt(chi[j][i-2] - chi[j][i-1]) * (2.0/3.0) * (xx[i-1] - xx[i-2]);
                    } else if (k == upper[j][i-1]) {
                        vx[j][i-1][k] = 0.0;
                    } else {
                        vx[j][i-1][k] = sqrt(2.0*(Uperp[j][k] - chi[j][i-1]));
                        mu[j][k] += (1.0/M_PI) * 0.5 * (vx[j][i-1][k] + vx[j][i-2][k]) * (xx[i-1] - xx[i-2]);
                    }
                    if (mu[j][k] != mu[j][k]) {
                        printf("BEFORE: mu[%d][%d] is NAN, kdrop[%d] = %d\n", j, k, j, kdrop[j]);
                        exit(-1);
                    }
                }
            } else if ((crossed_min[j] == 1) && (crossed_max[j] == 1) && (chi[j][i-1] < chiMax[j]) && (i-1 != imin[j])) {
                for (k = 0; k <= upperlimit[j]; k++) {
                    if ((chi[j][i-1] < Uperp[j][k]) && (chi[j][i-2] < Uperp[j][k])) {
                        vx[j][i-1][k] = sqrt(2.0*(Uperp[j][k] - chi[j][i-1]));
                        mu[j][k] += (1.0/M_PI) * 0.5 * (vx[j][i-1][k] + vx[j][i-2][k]) * (xx[i-1] - xx[i-2]);
                    } else if (Uperp[j][k] <= chi[j][i-1] && Uperp[j][k-1] > chi[j][i-1]) {
                        upper[j][i-1] = k;
                        ind = 0;
                        while (Uperp[j][k] < chi[j][i-2-ind]) ind++;
                        mu[j][k] += (sqrt(2.0)/M_PI) * (2.0/3.0) * (xx[i-1-ind] - xx[i-2-ind]) * pow(Uperp[j][k] - chi[j][i-2-ind], 1.5) / (chi[j][i-1-ind] - chi[j][i-2-ind]);
                    } else if (Uperp[j][k] <= chi[j][i-1] && Uperp[j][k] > chi[j][i-2]) {
                        mu[j][k] += (sqrt(2.0)/M_PI) * (2.0/3.0) * (xx[i-1] - xx[i-2]) * pow(Uperp[j][k] - chi[j][i-2], 1.5) / (chi[j][i-1] - chi[j][i-2]);
                    }
                    if (mu[j][k] != mu[j][k]) {
                        printf("mu[%d][%d] is NAN, kdrop[%d] = %d\n", j, k, j, kdrop[j]);
                        exit(-1);
                    }
                }
            } else if ((crossed_min[j] == 1) && (crossed_max[j] == 1) && (chi[j][i-1] > chiMax[j])) {
                xtop[j] = xx[i-2] + ((chiMax[j] - chi[j][i-2]) / (chi[j][i-1] - chi[j][i-2])) * (xx[i-1] - xx[i-2]);
                for (k = 0; k < upper[j][i-2]; k++) {
                    ind = 0;
                    while (Uperp[j][k] < chi[j][i-2-ind]) ind++;
                    mu[j][k] += (sqrt(2.0)/M_PI) * (2.0/3.0) * (xx[i-1-ind] - xx[i-2-ind]) * pow(Uperp[j][k] - chi[j][i-2-ind], 1.5) / (chi[j][i-1-ind] - chi[j][i-2-ind]);
                }
                crossed_max[j] = 0;
            }
            if (j != 0) {
                if (((chiMax[j-1] < chi[j-1][i-1] + TINY) && (chiMax[j] > chi[j][i-1] - TINY))) {
                    jmclosed[i-1] = j - 1;
                    if (i - 1 > icrit) {
                        jmopen[i-1] = j - 1;
                        if (DEBUG == 1) {
                            printf("i = %d, j = %d, icrit = %d\n", i-1, j-1, icrit);
                            intdvx_corr_delta = 0.0; intdvx_corr_delta_old = 0.0;
                            printf("jmopen[%d] = %d\n", i-1, jmopen[i-1]);
                        }
                    }
                }
            }
        }
    }

    upperlimit[sizexbar-1] = upperlimit[sizexbar-2] + 1;

    if (charge < 0 && fjmc_out != NULL) {
        fprintf(fjmc_out, "# x xbar_jmclosed xbar_jmopen\n");
        for (i = 0; i < size_finegrid; i++)
            fprintf(fjmc_out, "%f %f %f\n", xx[i], xbar[jmclosed[i]], xbar[jmopen[i]]);
    }

    maxj = sizexbar + 1;
    muopen       = calloc(maxj, sizeof(double));
    chiMopen     = calloc(maxj, sizeof(double));
    xbaropen     = calloc(maxj, sizeof(double));
    openorbitopen = calloc(maxj, sizeof(double));
    Ucritf       = calloc(maxj, sizeof(double));

    xbaropen[0]     = xbarcrit;
    chiMopen[0]     = chiMcrit;
    muopen[0]       = 0.0;
    openorbitopen[0] = 0.0;
    Ucritf[0]       = chiMcrit;
    printf("sizexbar = %d\n", sizexbar);

    for (j = 0; j < maxj - 1; j++) {
        mu[j][upperlimit[j]] = 0.0;
        chiMopen[j+1] = chiMax[j];
        k = 0;
        if (APPROXMUFORSMALLPHI == 1) {
            for (k = 0; k < upperlimit[j] + 1; k++) {
                if (fabs(phi[0]) < minphiformucalc) {
                    if (xbar[j] > 0.0) {
                        phibar   = lin_interp(xx, phi, xbar[j], size_finegrid, 1);
                        Omegaell = sqrt(1.0 + lin_interp(xx, phipp, xbar[j], size_finegrid, 1));
                        muell    = Uperp[j][upperlimit[j]-k] - phibar + 0.5 * pow(lin_interp(xx, phip, xbar[j], size_finegrid, 1) / Omegaell, 2.0);
                        muell   /= Omegaell;
                        musmall  = 0.5 * pow(xx[imin[j]] - xx[imin[j]-k], 2.0) * pow(chimpp[j], 0.5);
                        if (muapproxtype == 1) mu[j][upperlimit[j]-k] = muell;
                        else                   mu[j][upperlimit[j]-k] = musmall;
                    }
                }
            }
        }
        
        muopen[j+1]       = mu[j][0];
        xbaropen[j+1]     = xbar[j];
        Ucritf[j+1]       = chiMax[j] - mu[j][0];
        if (DEBUG == 1) printf("jmopen[%d] = %d\n", j, jmopen[j]);
        if (j == 0)
            openorbit[j] = (2.0*M_PI) * (((xbar[j] - xbarcrit) / (xbar[j+1] - xbarcrit)) * (mu[j+1][0] - mu[j][0]) / (xbar[j+1] - xbar[j]) + ((xbar[j+1] - xbar[j]) / (xbar[j+1] - xbarcrit)) * (mu[j][0] - 0.0) / (xbar[j] - xbarcrit));
        else if (j == maxj - 2)
            openorbit[j] = openorbit[j-1];
        else
            openorbit[j] = (2.0*M_PI) * (((xbar[j] - xbar[j-1]) / (xbar[j+1] - xbar[j-1])) * (mu[j+1][0] - mu[j][0]) / (xbar[j+1] - xbar[j]) + ((xbar[j+1] - xbar[j]) / (xbar[j+1] - xbar[j-1])) * (mu[j][0] - mu[j-1][0]) / (xbar[j] - xbar[j-1]));
        openorbitopen[j+1] = openorbit[j];
        openorbitantycal   = 2.0 * M_PI * xbar[j];
        if (DEBUG == 1) printf("%f %f %f %f %f\n", xbar[j], mu[j][0], Uperp[j][0], openorbit[j], openorbitantycal);
        if (TESTELL == 1) {
            printf("TESTELL entered\n");
            if (xbar[j] < xx[size_finegrid-1]) {
                printf("xbar condition fulfilled\n");
                for (k = 0; k < upperlimit[j]; k++) {
                    Omegaell = sqrt(1.0 + lin_interp(xx, phipp, xbar[j], size_finegrid, 1));
                    muell    = Uperp[j][k] - lin_interp(xx, phi, xbar[j], size_finegrid, 1) + 0.5 * pow(lin_interp(xx, phip, xbar[j], size_finegrid, 1) / Omegaell, 2.0);
                    muell   /= Omegaell;
                    if ((mu[j][k] < 2.5) && (mu[j][k] > 1.5))
                        printf("xbar = %f\tOmegaell = %f\tmu (actual, ellipmodel) = (%f, %f)\n", xbar[j], Omegaell, mu[j][k], muell);
                }
            }
        }
    }
    muopen[maxj-1]       = 999999.0;
    chiMopen[maxj-1]     = 999999.0;
    xbaropen[maxj-1]     = 100.0;
    openorbitopen[maxj-1] = 100.0;
    Ucritf[maxj-1]       = 0.0;

    if (DEBUG == 1) {
        printf("~~~~~The second element of FF is %f~~~~~\n", FF[0][1]);
        printf("~~~~~The second element of UU is %f~~~~~\n", UU[1]);
        printf("~~~~~The fourth element of mu is %f~~~~~\n", mumu[3]);
    }

    i = 0;
    clock_t int1 = clock();
    double inttime = (double)(int1 - begin) / CLOCKS_PER_SEC;
    if (DEBUG == 1) printf("in densfinorb_par: Array filling DONE: time is %f\n", inttime);

    if (fmu != NULL && charge < 0) {
        fprintf(fmu, "# j k xbar Uperp mu chiMax upperlimit\n");
        for (j = 0; j < sizexbar; j++) {
            for (k = 0; k <= upperlimit[j]; k++)
                fprintf(fmu, "%d %d %.15e %.15e %.15e %.15e %d\n",
                        j, k, xbar[j], Uperp[j][k], mu[j][k], chiMax[j], upperlimit[j]);
            fprintf(fmu, "\n");
        }
    }

    /* ================================================================
     * PHASE 2 — n_inf COMPUTATION  (identical to densfinorb)
     * ================================================================ */
    FILE *fout;
    if (Ti < 10.1) fout = fopen("OUTPUT/densfinorb_out.txt", "w");
    else           fout = fopen("TESTS/densfinorb_out.txt",  "w");
    if (fout == NULL) { printf("Cannot open densfinorb_out.txt"); exit(EXIT_FAILURE); }
    if (charge < 0) {
        fout_inref = fopen("OUTPUT/ne_inref.txt", "w");
        if (fout_inref == NULL) { printf("Cannot open ne_inref.txt\n"); exit(EXIT_FAILURE); }
    }

    ic = 0;
    while (x_grid[ic+1] < x_grid[size_phigrid-1] - limit_rho) ic++;
    size_xlim = ic;
    printf("size_xlim = %d/%d, limit_rho = %f\n", size_xlim, size_phigrid, limit_rho);
    deltax_inf = xx[size_xlim * zoomfactor] - xx[size_xlim * zoomfactor - 1];
    printf("deltax_inf = %f\n", deltax_inf);
    j_inf = (int)sqrt(2.0 * Ucap) / deltax_inf;
    printf("j_inf = %d, deltax_inf = %f\n", j_inf, deltax_inf);
    muinf  = (double*)calloc(j_inf, sizeof(double));
    vxinf  = (double**)calloc(j_inf, sizeof(double*));
    chiinf = (double*)calloc(j_inf, sizeof(double));
    for (j = j_inf - 1; j >= 0; j--) {
        vxinf[j] = (double*)calloc(j_inf, sizeof(double));
        chiinf[j] = 0.5 * deltax_inf * j * deltax_inf * j;
        for (k = j; k < j_inf; k++)
            vxinf[j][k] = sqrt(2.0 * (chiinf[k] - chiinf[j]));
        muinf[j] = chiinf[j];
    }

    intdxbar = 0.0;
    for (j = 0; j < j_inf; j++) {
        vxnew = 0.0;
        intdvxold = intdvx;
        intdvx = 0.0;
        intdU  = 0.0;
        for (k = j; k < j_inf; k++) {
            vxold    = vxnew;
            intdUold = intdU;
            intdU    = 0.0;
            munew    = muinf[k];
            vxnew    = vxinf[j][k];
            Ucrit    = lin_interp(muopen, Ucritf, munew, maxj, 790);
            Uperpnew = chiinf[k];
            sizeU    = (int)sqrt(2.0 * Ucap - 2.0 * Uperpnew) / dvz;
            reflected = (phi[0] > 0.0) ? 1 : 0;
            for (l = 0; l < sizeU; l++) {
                if (l != 0) {
                    Fold = F; Fold_ref = F;
                    vz = dvz * l;
                    U  = Uperpnew + 0.5 * vz * vz;
                    frac = 1.0;
                    F = bilin_interp(munew, U - munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                    if ((U - munew < Ucrit - numb) && (reflected == 1))
                        frac_reflected = 1.0;
                    else if ((U - munew > Ucrit - numb) && (reflected == 1)) {
                        reflected = 0;
                        if (Ucrit < Uperpnew - munew + numb) frac_reflected = 0.0;
                        else {
                            vzcrit    = sqrt(2.0 * (Ucrit + munew - Uperpnew));
                            Fold_ref  = bilin_interp(munew, Ucrit, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                            frac_reflected = (vzcrit - (vz - dvz)) / dvz;
                        }
                    } else frac_reflected = 0.0;
                    if (charge > 0.0) frac_reflected = 0.0;
                    intdU += 0.5 * frac * dvz * ((F + Fold) + frac_reflected * (F + Fold_ref));
                } else {
                    U = Uperpnew;
                    F = bilin_interp(munew, U - munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                    intdU += 0.0;
                }
            }
            intdUantycal = exp(-chiinf[k]) * (1.0 / (2.0 * M_PI));
            if ((j == 0) && (DEBUG == 1))
                printf("Analytical intdU is %f, numerical one is %f\n", intdUantycal, intdU);
            if (k != j) {
                dvx = vxnew - vxold;
                intdvx += 2.0 * 0.5 * dvx * (intdU + intdUold);
            }
            intdvxantycal = sqrt(1.0 / (2.0 * M_PI)) * exp(-0.5 * deltax_inf * j * deltax_inf * j);
        }
        if (DEBUG == 1) printf("intdvx is %f, analytical one is %f\n", intdvx, intdvxantycal);
        dxbar = deltax_inf;
        if (j != 0) intdxbar += 0.5 * (intdvx + intdvxold) * dxbar;
    }
    n_inf = 2.0 * intdxbar;
    printf("n_inf = %f\n", n_inf);
    if (DEBUG == 1) {
        printf("for charge = %f\n", charge);
        printf("n_inf = %f\n", n_inf);
        if ((n_inf != n_inf) || (n_inf < TINY)) { printf("n_inf = %f\n", n_inf); exit(-1); }
    }

    /* ================================================================
     * PHASE 3 — DENSITY INTEGRAL
     * ================================================================ */

    /* 3a. Diagnostic file setup */
    FILE *fupper = NULL;
    int upper_printed5 = 0, upper_printed6 = 0;
    if (charge < 0) fupper = fopen("OUTPUT/upper_diag.txt", "w");

    int i_49 = 0, i_51 = 0;
    {
        double diff49 = fabs(xx[0] - 4.5), diff51 = fabs(xx[0] - 5.5);
        for (int ii = 1; ii < size_finegrid; ii++) {
            if (fabs(xx[ii] - 4.5) < diff49) { diff49 = fabs(xx[ii] - 4.5); i_49 = ii; }
            if (fabs(xx[ii] - 5.5) < diff51) { diff51 = fabs(xx[ii] - 5.5); i_51 = ii; }
        }
    }
    FILE *f49 = fopen("OUTPUT/Uperp_min_x49.txt", "w");
    FILE *f51 = fopen("OUTPUT/Uperp_min_x51.txt", "w");
    (void)i_49; (void)i_51; /* written in sequential open-orbit loop below */

    int js_phi_imin = -1;
    if (!phi_monotone && phi_imin >= 0 && sizexbar > 0) {
        double best_xdiff = fabs(xbar[0] - xx[phi_imin]);
        js_phi_imin = 0;
        for (int js = 1; js < sizexbar; js++) {
            double xdiff = fabs(xbar[js] - xx[phi_imin]);
            if (xdiff < best_xdiff) { best_xdiff = xdiff; js_phi_imin = js; }
        }
    }

    /* Per-j storage for closed-orbit results */
    double *intdvx_j     = malloc(sizexbar * sizeof(double));
    double *intdvx_cd_j  = malloc(sizexbar * sizeof(double));
    double *intdvx_cm_j  = malloc(sizexbar * sizeof(double));
    double *intdvx_in_j  = malloc(sizexbar * sizeof(double));
    double *intdvx_ref_j = malloc(sizexbar * sizeof(double));
    if (!intdvx_j || !intdvx_cd_j || !intdvx_cm_j || !intdvx_in_j || !intdvx_ref_j) {
        fprintf(stderr, "densfinorb_par: out of memory\n");
        exit(EXIT_FAILURE);
    }

    stop = 0;
    ic   = 0;
    while (stop == 0) {
        i = ic * zoomfactor;

        memset(intdvx_j,     0, sizexbar * sizeof(double));
        memset(intdvx_cd_j,  0, sizexbar * sizeof(double));
        memset(intdvx_cm_j,  0, sizexbar * sizeof(double));
        memset(intdvx_in_j,  0, sizexbar * sizeof(double));
        memset(intdvx_ref_j, 0, sizexbar * sizeof(double));

        intdxbar = intdxbaropen = intdxbaropenflow = intdxbaropenener = 0.0;
        intdxbar_corr_delta = intdxbar_corr_chiM = 0.0;
        intdxbar_in = intdxbar_ref = 0.0;

        /* 3b. Open-orbit accumulation (sequential) */
        vxnew = 0.0;
        intdUopenold = 0.0; intdUopen = 0.0;
        intdUopenflowold = 0.0; intdUopenflow = 0.0;
        intdUopenenerold = 0.0; intdUopenener = 0.0;
        oorbintgrd = oorbintgrdold = 0.0;
        oorbintgrdflow = oorbintgrdflowold = 0.0;
        oorbintgrdener = oorbintgrdenerold = 0.0;

        for (j = 0; j < sizexbar; j++) {
            vxnew = 0.0;
            intdUopenold     = intdUopen;
            intdUopenflowold = intdUopenflow;
            intdUopenenerold = intdUopenener;
            intdUopen = intdUopenflow = intdUopenener = 0.0;

            if (j == jmopen[i]) {
                oorbintgrd = oorbintgrdold = 0.0;
                oorbintgrdflow = oorbintgrdflowold = 0.0;
                sizeU = (int)sqrt(2.0 * (Ucap - chiMax[j])) / dvzopen;
                if (j == 0) {
                    munew = 0.0; openorbitnew = 0.0;
                } else {
                    munew       = mu[j+1][0] + ((chiMax[j+1] - chi[j+1][i]) / (chiMax[j+1] - chi[j+1][i] + chi[j][i] - chiMax[j])) * (mu[j][0] - mu[j+1][0]);
                    openorbitnew = openorbit[j+1] + ((chiMax[j+1] - chi[j+1][i]) / (chiMax[j+1] - chi[j+1][i] + chi[j][i] - chiMax[j])) * (openorbit[j] - openorbit[j+1]);
                }
                for (l = 0; l < sizeU; l++) {
                    oorbintgrdold     = oorbintgrd;
                    oorbintgrdflowold = oorbintgrdflow;
                    vz = dvzopen * l;
                    if (j != 0) {
                        chinew  = chi[j+1][i] + ((chiMax[j+1] - chi[j+1][i]) / (chiMax[j+1] - chi[j+1][i] + chi[j][i] - chiMax[j])) * (chi[j][i] - chi[j+1][i]);
                        vx0open = 0.0;
                    } else {
                        chinew = chi[j][i];
                        vx0open = (chiMax[j] > chinew) ? sqrt(2.0 * (chiMax[j] - chinew)) : 0.0;
                    }
                    U = chinew + 0.5 * vz * vz;
                    if ((U > munew) && (U - 0.5*vz*vz + 0.5*(vz-dvzopen)*(vz-dvzopen) < munew)) {
                        frac    = (vz - sqrt(2.0 * (munew - chinew))) / dvzopen;
                        Fopen   = bilin_interp(munew, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                        oorbintgrdold = (sqrt(vx0open*vx0open + 2.0*alpha*sqrt(2.0*(munew-chinew))*openorbitnew) - vx0open) * Fopen;
                        oorbintgrdflowold = 0.5 * alpha * sqrt(2.0*(munew-chinew)) * openorbitnew * Fopen;
                        Fopen   = bilin_interp(munew, U - munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                        oorbintgrd     = (sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbitnew) - vx0open) * Fopen;
                        oorbintgrdflow = 0.5 * alpha * vz * openorbitnew * Fopen;
                    } else {
                        frac    = 1.0;
                        Fopen   = bilin_interp(munew, U - munew, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                        oorbintgrd     = (sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbitnew + TINY) - vx0open) * Fopen;
                        oorbintgrdflow = 0.5 * alpha * vz * openorbitnew * Fopen;
                    }
                    if (oorbintgrd != oorbintgrd) {
                        printf("vx0open = %f, openorbit[%d] = %f, imaginary oorbintgrd\n", vx0open, j, openorbit[j]);
                        exit(-1);
                    }
                    if (l != 0) {
                        intdUopen     += 0.5 * frac * dvzopen * (oorbintgrd     + oorbintgrdold);
                        intdUopenflow += 0.5 * frac * dvzopen * (oorbintgrdflow + oorbintgrdflowold);
                    }
                }
                if (j == 0) {
                    intdUopenold = 0.0;
                    dxbar        = 0.0;
                    intdxbaropen     += 0.5 * (intdUopen     + intdUopenold) * dxbar;
                    intdxbaropenflow += 0.5 * (intdUopenflow + intdUopenflowold) * dxbar;
                }
            } else if (j > jmopen[i]) {
                oorbintgrd = oorbintgrdold = 0.0;
                oorbintgrdflow = oorbintgrdflowold = 0.0;
                oorbintgrdener = oorbintgrdenerold = 0.0;
                sizeU = (int)sqrt(2.0 * (Ucap - chiMax[j])) / dvzopen;
                for (l = 0; l < sizeU; l++) {
                    oorbintgrdold     = oorbintgrd;
                    oorbintgrdflowold = oorbintgrdflow;
                    oorbintgrdenerold = oorbintgrdener;
                    vz      = dvzopen * l;
                    U       = chiMax[j] + 0.5 * vz * vz;
                    vx0open = sqrt(2.0 * (chiMax[j] - chi[j][i]));
                    if (vx0open != vx0open) {
                        printf("HERE imaginary vx0open, j = %d, i is %d, chi[j][i] = %f, chiMax[j] = %f\n", j, i, chi[j][i], chiMax[j]);
                        exit(-1);
                    }
                    if ((U >= mu[j][0]) && (U - 0.5*vz*vz + 0.5*(vz-dvzopen)*(vz-dvzopen) < mu[j][0])) {
                        frac     = (vz - sqrt(2.0 * (mu[j][0] - chiMax[j]))) / dvzopen;
                        Fopen    = bilin_interp(mu[j][0], 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                        oorbintgrdold     = sqrt(2.0 * alpha * sqrt(2.0*(mu[j][0]-chiMax[j])) * openorbit[j]) * Fopen;
                        oorbintgrdflowold = alpha * sqrt(2.0*(mu[j][0]-chiMax[j])) * openorbit[j] * Fopen;
                        oorbintgrdenerold = 0.125 * pow(2.0*alpha*sqrt(2.0*(mu[j][0]-chiMax[j]))*openorbit[j], 2.0) * Fopen + U * oorbintgrdflowold;
                        Fopen    = bilin_interp(mu[j][0], U - mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                        oorbintgrd     = sqrt(2.0 * alpha * vz * openorbit[j]) * Fopen;
                        oorbintgrdflow = alpha * vz * openorbit[j] * Fopen;
                        oorbintgrdener = 0.125 * pow(2.0*alpha*vz*openorbit[j], 2.0) * Fopen + U * oorbintgrdflow;
                    } else {
                        frac    = 1.0;
                        Fopen   = bilin_interp(mu[j][0], U - mu[j][0], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                        oorbintgrd     = (sqrt(vx0open*vx0open + 2.0*alpha*vz*openorbit[j]) - vx0open) * Fopen;
                        oorbintgrdflow = alpha * vz * openorbit[j] * Fopen;
                        oorbintgrdener = 0.125 * pow(2.0*alpha*vz*openorbit[j], 2.0) * Fopen + U * oorbintgrdflow;
                    }
                    if (oorbintgrd != oorbintgrd) {
                        printf("vx0open = %f, openorbit[%d] = %f, HERE imaginary oorbintgrd\n", vx0open, j, openorbit[j]);
                        exit(-1);
                    }
                    if (l != 0) {
                        intdUopen     += 0.5 * frac * dvzopen * (oorbintgrd     + oorbintgrdold);
                        intdUopenflow += 0.5 * frac * dvzopen * (oorbintgrdflow + oorbintgrdflowold);
                        intdUopenener += 0.5 * frac * dvzopen * (oorbintgrdener + oorbintgrdenerold);
                    }
                }
                if (intdUopen != intdUopen) { printf("HERE, j is %d\n", j); exit(-1); }
                dxbar = xbar[j] - xbar[j-1];
                if ((j == jmopen[i] + 1) && (j != 1))
                    dxbar = (xbar[j] - xbar[j-1]) * (chiMax[j] - chi[j][i]) / (chiMax[j] - chi[j][i] + chi[j-1][i] - chiMax[j-1]);
                intdxbaropen     += 0.5 * (intdUopen     + intdUopenold)     * dxbar;
                intdxbaropenflow += 0.5 * (intdUopenflow + intdUopenflowold) * dxbar;
                intdxbaropenener += 0.5 * (intdUopenener + intdUopenenerold) * dxbar;
                if (intdxbaropen != intdxbaropen) { printf("PROBLEM HERE, j is %d\n", j); exit(-1); }
            }
        } /* end sequential open-orbit j loop */

        /* 3c. Parallel closed-orbit k-loop */
        int jc0 = jmclosed[i];
        /* jc0 contributes intdvx = 0 */
        intdvx_j[jc0] = intdvx_cd_j[jc0] = intdvx_cm_j[jc0] = 0.0;
        intdvx_in_j[jc0] = intdvx_ref_j[jc0] = 0.0;

#pragma omp parallel for schedule(dynamic)
        for (j = jc0 + 1; j < sizexbar; j++) {
            /* All locals are thread-private */
            int    kk, ll, indd, sizeUU2, refl;
            double loc_dvx = 0.0;
            double loc_cd  = 0.0;
            double loc_cm  = 0.0;
            double loc_in  = 0.0;
            double loc_ref = 0.0;

            double lintdU = 0.0, lintdUold = 0.0;
            double lintdU_corr_delta = 0.0, lintdU_corr_delta_old = 0.0;
            double lintdU_corr_chiM  = 0.0, lintdU_corr_chiM_old  = 0.0;
            double lintdU_in = 0.0, lintdU_in_old = 0.0;
            double lintdU_ref = 0.0, lintdU_ref_old = 0.0;

            double lmunew = 0.0, lmuold = 0.0;
            double lUperpnew = 0.0, lUperpold = 0.0, lU_lb = 0.0;
            double lUcrit = 0.0, lUcritold = 0.0, lUcritp = 0.0, lUcritpold = 0.0;
            double ldmu_dUperp = 0.0, ldmu_dUperp_old = 0.0;
            double lU = 0.0, lvz = 0.0;
            double lF = 0.0, lFold = 0.0, lFold_ref = 0.0;
            double lfrac = 1.0, lfrac_reflected = 0.0, lvzcrit = 0.0;
            double lvxold = 0.0, lvxnew = 0.0, ldvx = 0.0;
            double lvx0open = 0.0;
            double lUperpnew_pre = 0.0, lUperp_lb_k = 0.0;
            double lmin_Uperp_j = Ucap, lmin_Uperp_bfx_j = Ucap, lUperp_lb_at_min_j = Ucap;
            double lintdUantycal = 0.0, lintdU_corr_delta2 = 0.0;

            for (kk = lowerlimit[j]; kk < upper[j][i] + 1; kk++) {
                lvxold        = lvxnew;
                lintdUold     = lintdU;
                lintdU        = 0.0;
                lintdU_corr_delta_old = lintdU_corr_delta;
                lintdU_corr_delta     = 0.0;
                lintdU_corr_chiM_old  = lintdU_corr_chiM;
                lintdU_corr_chiM      = 0.0;
                lintdU_in_old  = lintdU_in;
                lintdU_in      = 0.0;
                lintdU_ref_old = lintdU_ref;
                lintdU_ref     = 0.0;

                lmuold     = lmunew;
                lUcritold  = lUcrit;
                lUcritpold = lUcritp;
                lUperpold  = lUperpnew;
                ldmu_dUperp_old = ldmu_dUperp;

                if (kk == upper[j][i]) {
                    lUperpnew = chi[j][i];
                    if (lowerlimit[j] == upper[j][i])
                        lmunew = mu[j][kk];
                    else
                        lmunew = ((chi[j][i] - Uperp[j][kk]) * mu[j][kk-1] + (Uperp[j][kk-1] - chi[j][i]) * mu[j][kk]) / (Uperp[j][kk-1] - Uperp[j][kk]);
                    lvxnew = 0.0;
                    ldvx   = lvxold - lvxnew;
                } else {
                    lUperpnew = Uperp[j][kk];
                    lvxnew    = vx[j][i][kk];
                    ldvx      = lvxold - lvxnew;
                    lmunew    = mu[j][kk];
                }

                lUcrit = lin_interp(muopen, Ucritf, lmunew, maxj, 1879);

                lUperpnew_pre = lUperpnew;
                lUperp_lb_k   = lUperpnew;
                lU_lb         = lUperpnew;
                if (!phi_monotone) {
                    double lUperp_lb = lUperpnew;
                    for (int js = j + 1; js <= (js_phi_imin + sizexbar) / 2; js++) {
                        int k_found = -1;
                        double Uperps = uperp_from_mu2(js, lmunew, mu, Uperp, lowerlimit[js], upper[js][i], upperlimit[js], j, kk, &k_found, phi[phi_imin]);
                        if (Uperps > lUperp_lb) lUperp_lb = Uperps;
                    }
                    if (lUperp_lb > lUperpnew) lU_lb = lUperp_lb;
                    lUperp_lb_k = lUperp_lb;
                }
                if (lUperpnew_pre < lmin_Uperp_bfx_j) {
                    lmin_Uperp_bfx_j   = lUperpnew_pre;
                    lUperp_lb_at_min_j = lUperp_lb_k;
                }
                if (lUperpnew < lmin_Uperp_j) lmin_Uperp_j = lUperpnew;

                sizeUU2 = (int)sqrt(2.0 * (Ucap - lU_lb)) / dvz;
                refl = 1;
                for (ll = 0; ll < sizeUU2; ll++) {
                    if (ll != 0) {
                        lFold = lF; lFold_ref = lF;
                        lU  = lU_lb + 0.5 * (dvz*ll) * (dvz*ll);
                        lvz = dvz * ll;
                        if ((lU > lmunew) && (lU - 0.5*lvz*lvz + 0.5*(lvz-dvz)*(lvz-dvz) < lmunew)) {
                            lfrac  = (lvz - sqrt(2.0 * (lmunew - lU_lb + TINY))) / dvz;
                            lFold  = bilin_interp(lmunew, 0.0, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                            lF     = bilin_interp(lmunew, lU - lmunew, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                        } else if (lU > lmunew) {
                            lfrac = 1.0;
                            lF    = bilin_interp(lmunew, lU - lmunew, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                        } else {
                            lfrac = 1.0;
                            lF    = 0.0;
                        }
                        if ((lU - lmunew < lUcrit - numb) && (refl == 1))
                            lfrac_reflected = 1.0;
                        else if ((lU - lmunew > lUcrit - numb) && (refl == 1)) {
                            refl = 0;
                            if (lUcrit < lU_lb - lmunew + numb) lfrac_reflected = 0.0;
                            else {
                                lvzcrit    = sqrt(2.0 * (lUcrit + lmunew - lU_lb));
                                lFold_ref  = bilin_interp(lmunew, lUcrit, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                                lfrac_reflected = (lvzcrit - (lvz - dvz)) / dvz;
                            }
                        } else {
                            lfrac_reflected = 0.0;
                        }
                        lintdU += 0.5 * lfrac * dvz * ((lF + lFold) + lfrac_reflected * (lF + lFold_ref));
                        if (lintdU != lintdU) {
                            printf("intdU is NAN, j=%d, i=%d\n", j, i);
                            exit(-1);
                        }
                        if (charge < 0) {
                            lintdU_in  += 0.5 * lfrac * dvz * (lF + lFold);
                            lintdU_ref += 0.5 * lfrac * dvz * lfrac_reflected * (lF + lFold_ref);
                        }
                        if (kk == lowerlimit[j] && charge < 0)
                            lintdU_corr_chiM = lintdU;
                    } else {
                        lU = lU_lb;
                        lF = bilin_interp(lmunew, lU - lmunew, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                        lintdU += 0.0;
                        if (kk == lowerlimit[j] && charge < 0)
                            lintdU_corr_chiM = lintdU;
                    }
                }
                lintdUantycal = exp(-lUperpnew) * (1.0 / (2.0 * M_PI));
                if (lUcrit + lmunew > lUperpnew)
                    lintdU_corr_delta2 = bilin_interp(lmunew, lUcrit, FF, mumu, UU, sizemumu, sizeUU, -1, -1);
                else
                    lintdU_corr_delta2 = 0.0;
                lintdU_corr_delta = lintdU_corr_delta2;

                if (kk == lowerlimit[j]) {
                    loc_dvx += 0.0;
                    loc_cd  += 0.0;
                    if (charge < 0) {
                        lvx0open = sqrt(2.0 * (chiMax[j] - chi[j][i]));
                        loc_cm   = 2.0 * lintdU_corr_chiM / lvx0open;
                    }
                } else {
                    double ldmu     = lmunew - lmuold;
                    double ldmu_tol = fmax(TINY, 1e-14 * fmax(1.0, fmax(fabs(lmunew), fabs(lmuold))));
                    double ldUperp     = lUperpnew - lUperpold;
                    double ldUperp_tol = fmax(TINY, 1e-14 * fmax(1.0, fmax(fabs(lUperpnew), fabs(lUperpold))));

                    loc_dvx += 2.0 * 0.5 * ldvx * (lintdU + lintdUold);
                    if (charge < 0) {
                        loc_in  += 2.0 * 0.5 * ldvx * (lintdU_in  + lintdU_in_old);
                        loc_ref += 2.0 * 0.5 * ldvx * (lintdU_ref + lintdU_ref_old);
                    }
                    if (fabs(ldmu) <= ldmu_tol) lUcritp = lUcritpold;
                    else lUcritp = (lUcrit - lUcritold) / ldmu;

                    if (fabs(ldUperp) <= ldUperp_tol) {
                        if (fabs(ldmu) <= ldmu_tol) ldmu_dUperp = 0.0;
                        else ldmu_dUperp = ldmu_dUperp_old;
                    } else {
                        ldmu_dUperp = ldmu / ldUperp;
                        if (fabs(ldmu_dUperp) > 10000) {
                            printf("At x = %f, xbar = %f, Uperp = %.6e, Uperpold = %.6e, munew = %.6e, muold = %.6e, munew - muold = %.6e, Uperpnew - Uperpold is %.6e (tol %.6e), clamping, k - upperlimit[j] = %d\n",
                                   xx[i], xbar[j], lUperpnew, lUperpold, lmunew, lmuold, ldmu, ldUperp, ldUperp_tol, kk - upperlimit[j]);
                            ldmu_dUperp = ldmu_dUperp_old;
                        }
                    }
                    loc_cd += 2.0 * 0.5 * ldvx * (lintdU_corr_delta * (1.0 + lUcritp * ldmu_dUperp) + lintdU_corr_delta_old * (1.0 + lUcritpold * ldmu_dUperp_old));
                }
                if (loc_dvx != loc_dvx) {
                    printf("intdvx is NAN, j=%d, i=%d\n", j, i);
                    exit(-1);
                }
                (void)lintdUantycal;
                (void)lmin_Uperp_j; (void)lmin_Uperp_bfx_j; (void)lUperp_lb_at_min_j;
            } /* end k loop */

            intdvx_j[j]     = loc_dvx;
            intdvx_cd_j[j]  = loc_cd;
            intdvx_cm_j[j]  = loc_cm;
            intdvx_in_j[j]  = loc_in;
            intdvx_ref_j[j] = loc_ref;
        } /* end #pragma omp parallel for */

        /* 3d. Sequential trapezoidal reduction over j */
        for (j = jc0 + 1; j < sizexbar; j++) {
            double dxbj;
            if (j == jc0 + 1 && i != 0) {
                if (i == imax[j])
                    dxbj = xbar[j] - xbar[j-1];
                else
                    dxbj = (xbar[j] - xbar[j-1]) * (chiMax[j] - chi[j][i]) / (chiMax[j] - chi[j][i] + chi[j-1][i] - chiMax[j-1]);
            } else {
                dxbj = xbar[j] - xbar[j-1];
            }
            intdxbar            += 0.5 * (intdvx_j[j]     + intdvx_j[j-1])     * dxbj;
            intdxbar_corr_delta += 0.5 * (intdvx_cd_j[j]  + intdvx_cd_j[j-1])  * dxbj;
            if (charge < 0) {
                intdxbar_corr_chiM += 0.5 * (intdvx_cm_j[j]  + intdvx_cm_j[j-1])  * dxbj;
                intdxbar_in        += 0.5 * (intdvx_in_j[j]  + intdvx_in_j[j-1])  * dxbj;
                intdxbar_ref       += 0.5 * (intdvx_ref_j[j] + intdvx_ref_j[j-1]) * dxbj;
            }
            if (intdxbar != intdxbar) {
                printf("intdxbar is NAN, j=%d, i=%d\n", j, i);
                exit(-1);
            }
        }

        /* 3e. Diagnostic upper_diag, n_grid write, stop logic */
        if (fupper != NULL) {
            if (!upper_printed5 && xx[i] >= 4.0) {
                upper_printed5 = 1;
                FILE *f5 = fopen("OUTPUT/upper_diag_x5.txt", "w");
                if (f5 != NULL) {
                    fprintf(f5, "# x=%.6f ic=%d (target 5.0)\n", xx[i], ic);
                    fprintf(f5, "# j xbar upper Uperp_inner chi mu\n");
                    for (j = 0; j < sizexbar; j++)
                        fprintf(f5, "%d %.6f %d %.6f %.6f %.6e\n", j, xbar[j], upper[j][i], Uperp[j][upper[j][i]], chi[j][i], mu[j][upper[j][i]]);
                    fclose(f5);
                }
            }
            if (!upper_printed6 && xx[i] >= 6.0) {
                upper_printed6 = 1;
                FILE *f6 = fopen("OUTPUT/upper_diag_x6.txt", "w");
                if (f6 != NULL) {
                    fprintf(f6, "# x=%.6f ic=%d (target 6.0)\n", xx[i], ic);
                    fprintf(f6, "# j xbar upper Uperp_inner chi mu\n");
                    for (j = 0; j < sizexbar; j++)
                        fprintf(f6, "%d %.6f %d %.6f %.6f %.6e\n", j, xbar[j], upper[j][i], Uperp[j][upper[j][i]], chi[j][i], mu[j][upper[j][i]]);
                    fclose(f6);
                }
            }
        }

        n_grid[ic] = intdxbar + intdxbaropen;
        n_grid_corr_delta[ic] = intdxbar_corr_delta;
        if (charge < 0) n_grid_corr_chiM[ic] = intdxbar_corr_chiM;

        if (ic == 0) {
            flux0  = intdxbaropenflow / (n_inf * alpha);
            Qflux0 = intdxbaropenener / (n_inf * alpha);
            if (charge < 0.0) *flux = flux0;
            printf("flow velocity at x=0 = %f\n", flux0 / n_grid[0]);
            printf("flux evaluated at x=0 is %f\theat flux = %f\n", flux0, Qflux0);
        }
        if ((ic != size_xlim) && (1.0 - n_grid[ic]/n_inf <= margin) && (ic > 0 && phi_grid[ic] > phi_grid[ic-1]) && phi_grid[ic] < 0.0) {
            printf("ic = %d/%d, size_xlim = %d, margin = %f, condition (< margin?) = %f\n", ic, size_phigrid, size_xlim, margin, 1.0 - n_grid[ic]/n_inf);
            stop = 1;
            *size_ngrid = ic;
            printf("stopping density evaluation at x = %f, density = %f*n_inf\n", x_grid[ic], n_grid[ic]/n_inf);
        }
        if (DEBUG == 1) {
            printf("%f, %f, %f is TOTAL, CLOSED and OPEN orbit density at position index %d, position %f, potential = %f, n_inf %f, ic = %d, size_xlim = %d\n",
                   n_grid[ic], intdxbar, intdxbaropen, ic, xx[i], phi_grid[ic], n_inf, ic, size_xlim);
        }
        if (ic == size_xlim - 1) {
            printf("ic = %d/%d\n", ic, size_phigrid);
            stop = 1;
            *size_ngrid = ic;
            printf("In densfinorb_par.c: reached maximum distance from wall\n");
        }
        if (ic == size_xlim) {
            n_inf = n_grid[ic];
            ic = 0;
        } else {
            fprintf(fout, "%f %f %f %f\n", xx[i], n_grid[ic]/n_inf, intdxbar/n_inf, intdxbaropen/n_inf);
            if (fout_inref != NULL)
                fprintf(fout_inref, "%f %f %f %f\n", xx[i], n_grid[ic]/n_inf, intdxbar_in/n_inf, intdxbar_ref/n_inf);
            ic += 1;
        }
    } /* end while (stop == 0) */

    /* ================================================================
     * POST-PROCESSING
     * ================================================================ */
    printf("NINF = %f\n", n_inf);
    printf("in densfinorb_par: charge = %f, dndphi = %f\n", charge, (n_grid[*size_ngrid-1] - n_grid[*size_ngrid-2]) / (phi_grid[*size_ngrid-1] - phi_grid[*size_ngrid-2]));
    fclose(fout);
    if (fout_inref != NULL) fclose(fout_inref);
    if (fupper != NULL) fclose(fupper);
    if (f49 != NULL) fclose(f49);
    if (f51 != NULL) fclose(f51);
    if (stop == 0) { printf("ERROR: the density never reached stopdens*n_inf\n"); exit(-1); }

    for (ic = 0; ic < size_phigrid; ic++) {
        if (ic < *size_ngrid) {
            n_grid[ic]            /= n_inf;
            n_grid_corr_delta[ic] /= n_inf;
            if (charge < 0) n_grid_corr_chiM[ic] /= n_inf;
        } else {
            n_grid[ic] = 0.0;
        }
    }

    for (j = 0; j < maxj/2; j++) {
        vy_op[j]     = xbaropen[j];
        mu_op[j]     = muopen[j];
        chiMax_op[j] = chiMopen[j];
        dmudvy_op[j] = openorbitopen[j];
        if (DEBUG == 1)
            printf("%d/%d %f %f %f %f\n", j, sizemumu, mu_op[j], vy_op[j], chiMax_op[j], dmudvy_op[j]);
    }
    *size_op = maxj / 2;

    if (filellip != NULL) fclose(filellip);

    free(chiMopen);
    free(openorbitopen);
    free(xbaropen);
    free(Ucritf);
    free(muopen);
    free(xx);
    free(jmclosed);
    free(jmopen);
    free(xifunction);
    free(xbar);
    free(chimin);
    free(chiMax);
    free(chimpp);
    free(crossed_min);
    free(crossed_max);
    free(openorbit);
    free(upperlimit);
    free(lowerlimit);
    free(xtop);
    free(imax);
    free(imin);

    for (w = 0; w < sizexbar; w++) {
        free(chi[w]);
        free(Uperp[w]);
        free(mu[w]);
        free(upper[w]);
        for (s = 0; s < size_finegrid; s++) free(vx[w][s]);
        free(vx[w]);
    }
    free(chi);
    free(Uperp);
    free(mu);
    free(upper);
    free(vx);

    free(phi);
    free(gg);
    free(ff);
    free(phip);
    free(kdrop);

    for (j = 0; j < j_inf; j++) free(vxinf[j]);
    free(vxinf);
    free(chiinf);
    free(muinf);

    free(intdvx_j);
    free(intdvx_cd_j);
    free(intdvx_cm_j);
    free(intdvx_in_j);
    free(intdvx_ref_j);

    clock_t end = clock();
    double jobtime = (double)(end - begin) / CLOCKS_PER_SEC;
    double wt_end  = omp_get_wtime();
    printf("in densfinorb_par: module ran in %.2f s CPU time, %.2f s wall time (%d thread(s))\n",
           jobtime, wt_end - wt_begin, omp_get_max_threads());
    return;
}
