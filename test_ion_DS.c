/*
 * test_ion_DS.c — convergence test for ion DS density (densionDS2)
 *
 * Usage: ./test_ion_DS <phi_MP_file> <phi_DS_file> <alpha_deg> <TiovTe> <margin_mp>
 *                      [zoom1 zoom2 ...]
 *
 * phi_MP_file : two-column text file (x  phi) for the magnetic presheath.
 *               x=0 is the DSE (inner wall of the MP); x increases toward
 *               upstream (phi → 0).  phi[0] is the DSE potential and is
 *               used as phi0 in densionDS2.
 * phi_DS_file : two-column text file (x  phi) for the Debye sheath.
 *               x=0 is the DS wall; x increases toward the DSE (phi → 0).
 * alpha_deg   : magnetic-field angle alpha in degrees
 * TiovTe      : ion-to-electron temperature ratio
 * margin_mp   : stopping margin for the MP densfinorb call (e.g. 0.04).
 *               Use a negative value to always run to the end of the grid.
 * zoomN       : zoom factors for the MP spatial grid (default: 1 2 4 8)
 *
 * For each zoom z:
 *   1. Runs densfinorb (charge=+1) on phi_MP to get open-orbit quantities
 *      (vy, mu_op, chiM, twopidmudvy, size_op).
 *   2. Runs densionDS2 on phi_DS with phi0 = phi_MP[0].
 *   3. Writes OUTPUT/test_ion_DS_ni_zoom<z>.txt  (columns: x  phi  ni)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include "mps.h"

#define MAXMU_I   8.0
#define MAXVPAR_I 5.0
#define DMU_I     0.05
#define DVPAR_I   0.05
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void read_phi_file(const char *path, double **px, double **pphi, int *pn)
{
    int cap = 256, n = 0;
    double *x   = malloc(cap * sizeof(double));
    double *phi = malloc(cap * sizeof(double));
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", path); exit(1); }
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == ' ') continue;
        double xv, pv;
        if (sscanf(line, "%lf %lf", &xv, &pv) != 2) continue;
        if (n >= cap) {
            cap *= 2;
            x   = realloc(x,   cap * sizeof(double));
            phi = realloc(phi, cap * sizeof(double));
        }
        x[n] = xv; phi[n] = pv; n++;
    }
    fclose(fp);
    *px = x; *pphi = phi; *pn = n;
}

/* Append a constant-phi tail so densfinorb's internal limit_rho=8 cutoff
 * does not truncate the original n points. */
static void make_ext(const double *x, const double *phi, int n,
                     double **x_out, double **phi_out, int *n_out)
{
    double dx = x[n-1] - x[n-2];
    int n_tail = (int)ceil(9.0 / dx) + 2;
    int ne = n + n_tail;
    double *xe   = malloc(ne * sizeof(double));
    double *phie = malloc(ne * sizeof(double));
    memcpy(xe,   x,   n * sizeof(double));
    memcpy(phie, phi, n * sizeof(double));
    for (int i = 0; i < n_tail; i++) {
        xe[n + i]   = x[n-1]   + (i + 1) * dx;
        phie[n + i] = phi[n-1];
    }
    *x_out = xe; *phi_out = phie; *n_out = ne;
}

int main(int argc, char *argv[])
{
    if (argc < 6) {
        fprintf(stderr,
            "Usage: %s <phi_MP_file> <phi_DS_file> <alpha_deg> <TiovTe>"
            " <margin_mp> [zoom1 zoom2 ...]\n", argv[0]);
        return 1;
    }
    const char *phi_MP_file = argv[1];
    const char *phi_DS_file = argv[2];
    double alpha     = atof(argv[3]) * M_PI / 180.0;
    double TiovTe    = atof(argv[4]);
    double margin_mp = atof(argv[5]);

    int nzooms = (argc > 6) ? argc - 6 : 4;
    int *zooms = malloc(nzooms * sizeof(int));
    if (argc > 6) {
        for (int z = 0; z < nzooms; z++) zooms[z] = atoi(argv[6 + z]);
    } else {
        int def[] = {1, 2, 4, 8};
        for (int z = 0; z < 4; z++) zooms[z] = def[z];
    }

    double *x_MP, *phi_MP, *x_DS, *phi_DS;
    char ds_base[256];
    {
        const char *s = strrchr(phi_DS_file, '/');
        strncpy(ds_base, s ? s + 1 : phi_DS_file, sizeof(ds_base) - 1);
        ds_base[sizeof(ds_base) - 1] = '\0';
        char *dot = strrchr(ds_base, '.');
        if (dot) *dot = '\0';
    }

    int n_MP, n_DS;
    read_phi_file(phi_MP_file, &x_MP, &phi_MP, &n_MP);
    read_phi_file(phi_DS_file, &x_DS, &phi_DS, &n_DS);
    printf("MP: %d points  x=[%.4f, %.4f]  phi=[%.6f, %.6f]\n",
           n_MP, x_MP[0], x_MP[n_MP-1], phi_MP[0], phi_MP[n_MP-1]);
    printf("DS: %d points  x=[%.4f, %.4f]  phi=[%.6f, %.6f]\n",
           n_DS, x_DS[0], x_DS[n_DS-1], phi_DS[0], phi_DS[n_DS-1]);

    double phi0 = phi_MP[0];
    printf("phi0 (DSE potential) = %.6f\n", phi0);

    double *x_MP_ext, *phi_MP_ext;
    int n_MP_ext;
    make_ext(x_MP, phi_MP, n_MP, &x_MP_ext, &phi_MP_ext, &n_MP_ext);

    /* Ion distribution via Figen2 (Chodura-condition shifted distribution) */
    int size_mu_i   = (int)(MAXMU_I  / DMU_I)   + 1;
    int size_vpar_i = (int)(MAXVPAR_I / DVPAR_I) + 1;
    double *mu_i    = malloc(size_mu_i   * sizeof(double));
    double *U_i     = malloc(size_vpar_i * sizeof(double));
    double **dist_i = malloc(size_mu_i   * sizeof(double *));
    for (int i = 0; i < size_mu_i; i++)
        dist_i[i] = malloc(size_vpar_i * sizeof(double));

    /* Figen2 expects species-indexed arrays; wrap single species here */
    int    sizevperp_arr[1] = { size_mu_i };
    int    sizevpar_arr[1]  = { size_vpar_i };
    double nioverne_arr[1]  = { 1.0 };
    double mioverme_arr[1]  = { 1.0 };
    double TioverTe_arr[1]  = { TiovTe };
    Figen2(&dist_i, &U_i, &mu_i, 1,
           nioverne_arr, mioverme_arr, TioverTe_arr,
           sizevpar_arr, sizevperp_arr,
           DVPAR_I, DMU_I);

    int max_op = n_MP_ext + 64;
    double *ni_MP      = calloc(n_MP_ext, sizeof(double));
    double *corr_d     = calloc(n_MP_ext, sizeof(double));
    double *corr_c     = calloc(n_MP_ext, sizeof(double));
    double *vy_op      = malloc(max_op * sizeof(double));
    double *mu_op      = malloc(max_op * sizeof(double));
    double *chiM       = malloc(max_op * sizeof(double));
    double *dmudvy     = malloc(max_op * sizeof(double));

    double *ni_DS      = calloc(n_DS, sizeof(double));
    double *ni_DS_corr = calloc(n_DS, sizeof(double));
    double *ni_DS_refl = calloc(n_DS, sizeof(double));

    mkdir("OUTPUT", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for (int z = 0; z < nzooms; z++) {
        int zoom = zooms[z];
        int size_ni_MP = 0, size_op = 0;
        double flux = 0.0, garbage = 0.0, Bohm = 0.0;

        memset(ni_MP,      0, n_MP_ext * sizeof(double));
        memset(corr_d,     0, n_MP_ext * sizeof(double));
        memset(corr_c,     0, n_MP_ext * sizeof(double));
        memset(ni_DS,      0, n_DS * sizeof(double));
        memset(ni_DS_corr, 0, n_DS * sizeof(double));
        memset(ni_DS_refl, 0, n_DS * sizeof(double));

        printf("\n=== zoom = %d ===\n", zoom);

        densfinorb(TiovTe, 1.0, alpha, n_MP_ext, &size_ni_MP,
                   ni_MP, corr_d, corr_c, x_MP_ext, phi_MP_ext, +1.0,
                   dist_i, mu_i, U_i, size_mu_i, size_vpar_i,
                   0.0, &flux, &garbage, zoom,
                   margin_mp, -999.9,
                   vy_op, mu_op, chiM, dmudvy, &size_op,
                   NULL, NULL);
        printf("MP densfinorb: size_ni_MP=%d  size_op=%d  flux=%.6f\n",
               size_ni_MP, size_op, flux);

        densionDS2(alpha, TiovTe, &Bohm,
                   ni_DS, phi_DS, phi0,
                   dist_i, mu_i, U_i,
                   vy_op, mu_op, chiM, dmudvy,
                   n_DS, size_mu_i, size_vpar_i, size_op,
                   ni_DS_corr, ni_DS_refl);
        printf("densionDS2: Bohm=%.6f\n", Bohm);

        char outpath[256];
        snprintf(outpath, sizeof(outpath),
                 "OUTPUT/test_ion_DS_ni_%s_zoom%d.txt", ds_base, zoom);
        FILE *fo = fopen(outpath, "w");
        if (!fo) { fprintf(stderr, "Cannot open %s\n", outpath); continue; }
        fprintf(fo, "# x phi ni ni_corr, ni_ref\n");
        for (int i = 0; i < n_DS; i++)
            fprintf(fo, "%.12e %.12e %.12e %.12e %.12e\n",
                    x_DS[i], phi_DS[i], ni_DS[i], ni_DS_corr[i], ni_DS_refl[i]);
        fclose(fo);
        printf("Wrote %s (%d points)\n", outpath, n_DS);
    }

    for (int i = 0; i < size_mu_i; i++) free(dist_i[i]);
    free(dist_i); free(mu_i); free(U_i);
    free(x_MP); free(phi_MP); free(x_DS); free(phi_DS);
    free(x_MP_ext); free(phi_MP_ext);
    free(ni_MP); free(corr_d); free(corr_c);
    free(vy_op); free(mu_op); free(chiM); free(dmudvy);
    free(ni_DS); free(ni_DS_corr); free(ni_DS_refl);
    free(zooms);
    return 0;
}
