/*
 * test_DS.c — standalone convergence test for densfinorb (DS electrons)
 *
 * Usage: ./test_DS <phi_file> <alpha> [zoom1 zoom2 ...]
 *
 * phi_file : two-column text file (x  phi), ordered from wall (x=0, most
 *            negative phi) to DS entrance (largest x, phi≈0). Lines starting
 *            with '#' are skipped.
 * alpha    : magnetic-field angle parameter (same as in GYRAZE)
 * zoomN    : integer zoom factors to test (default: 1 2 4 8)
 *
 * For each zoom factor z the internal fine grid has spacing deltax/z.
 * Output: OUTPUT/test_DS_ne_zoom<z>.txt  (columns: x  phi  ne)
 *
 * Distribution function: Maxwellian, F(mu, E_par) = exp(-(mu+E_par))/(2*pi)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <omp.h>
#include "mps.h"

#define MAXMU    10.0
#define MAXVPAR  10.0
#define DMU      0.1
#define DVPAR    0.05
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

int main(int argc, char *argv[])
{
    if (argc < 5) {
        fprintf(stderr,
            "Usage: %s <phi_DS_file> <phi_MP_file> <alpha> <margin> [--parallel] [zoom1 zoom2 ...]\n"
            "  phi_DS_file: DS potential (x=0 wall, phi->0 at DSE)\n"
            "  phi_MP_file: MP potential (x=0 DSE, phi->0 upstream); phi_MP[0] sets\n"
            "               the presheath energy shift for the electron distribution.\n"
            "  margin: stopping threshold (same as MARGIN_DS in GYRAZE).\n"
            "          Use a negative value (e.g. -1.0) to disable the\n"
            "          density-overshoot stop and always run to the end\n"
            "          of the phi grid.\n"
            "  --parallel: use OpenMP parallel densfinorb_par instead of serial densfinorb\n",
            argv[0]);
        return 1;
    }
    const char *phi_DS_file = argv[1];
    const char *phi_MP_file = argv[2];
    double alpha  = atof(argv[3])*M_PI/180;
    double margin = atof(argv[4]);

    int use_parallel = 0;
    int zoom_argc = 0;
    char **zoom_argv = malloc(argc * sizeof(char *));
    for (int i = 5; i < argc; i++) {
        if (strcmp(argv[i], "--parallel") == 0) use_parallel = 1;
        else zoom_argv[zoom_argc++] = argv[i];
    }
    printf("use_parallel = %d\n", use_parallel);
#define DENSFINORB(...) (use_parallel ? densfinorb_par(__VA_ARGS__) : densfinorb(__VA_ARGS__))

    int nzooms = (zoom_argc > 0) ? zoom_argc : 4;
    int *zooms = malloc(nzooms * sizeof(int));
    if (zoom_argc > 0) {
        for (int z = 0; z < nzooms; z++) zooms[z] = atoi(zoom_argv[z]);
    } else {
        int def[] = {1, 2, 4, 8};
        for (int z = 0; z < 4; z++) zooms[z] = def[z];
    }
    free(zoom_argv);

    char ds_base[256];
    {
        const char *s = strrchr(phi_DS_file, '/');
        strncpy(ds_base, s ? s + 1 : phi_DS_file, sizeof(ds_base) - 1);
        ds_base[sizeof(ds_base) - 1] = '\0';
        char *dot = strrchr(ds_base, '.');
        if (dot) *dot = '\0';
    }

    double *x_DS, *phi_DS;
    int n_phi;
    read_phi_file(phi_DS_file, &x_DS, &phi_DS, &n_phi);
    printf("Read %d points from %s\n", n_phi, phi_DS_file);
    printf("  x   range: [%.6f, %.6f]\n", x_DS[0],   x_DS[n_phi-1]);
    printf("  phi range: [%.6f, %.6f]\n", phi_DS[0], phi_DS[n_phi-1]);

    double *x_MP, *phi_MP;
    int n_MP;
    read_phi_file(phi_MP_file, &x_MP, &phi_MP, &n_MP);
    double phi_DSE = phi_MP[0];
    printf("phi_DSE = phi_MP[0] = %.6f  (presheath energy shift)\n", phi_DSE);
    free(x_MP); free(phi_MP);

    /* densfinorb stops limit_rho=8 units before the end of the grid.
     * Append a constant-phi tail so the full n_phi points get computed. */
    double dx_last = x_DS[n_phi-1] - x_DS[n_phi-2];
    int n_tail   = (int)ceil(9.0 / dx_last) + 2;
    int n_ext    = n_phi + n_tail;
    double *x_ext   = malloc(n_ext * sizeof(double));
    double *phi_ext = malloc(n_ext * sizeof(double));
    memcpy(x_ext,   x_DS,   n_phi * sizeof(double));
    memcpy(phi_ext, phi_DS, n_phi * sizeof(double));
    for (int i = 0; i < n_tail; i++) {
        x_ext[n_phi + i]   = x_DS[n_phi-1]  + (i + 1) * dx_last;
        phi_ext[n_phi + i] = phi_DS[n_phi-1];
    }

    /* Build Maxwellian in MP parallel-speed coordinates (dist_DK), then map
     * onto the DS energy grid using the presheath energy shift phi_DSE:
     *   vpar_MP = sqrt(2*U_DS - 2*phi_DSE)   (energy conservation, phi_DSE < 0)
     * This follows the same transformation as GYRAZE lines 2403-2408.
     * For a Maxwellian dist_DK the result equals exp(phi_DSE)*Maxwellian(U_DS),
     * a constant factor that cancels in densfinorb's normalization; the
     * transformation matters when dist_DK is non-Maxwellian (presheath-depleted). */
    int size_mu   = (int)(MAXMU  / DMU)   + 1;
    int size_vpar = (int)(MAXVPAR / DVPAR) + 1;
    double *mu_e    = malloc(size_mu   * sizeof(double));
    double *vpar_MP = malloc(size_vpar * sizeof(double));
    double *U_e_DS  = malloc(size_vpar * sizeof(double));
    double **dist_DK = malloc(size_mu  * sizeof(double *));
    double **dist    = malloc(size_mu  * sizeof(double *));
    for (int i = 0; i < size_mu; i++) {
        mu_e[i] = i * DMU;
        dist_DK[i] = malloc(size_vpar * sizeof(double));
        dist[i]    = malloc(size_vpar * sizeof(double));
        for (int j = 0; j < size_vpar; j++) {
            if (i == 0) {
                vpar_MP[j] = j * DVPAR;
                U_e_DS[j]  = 0.5 * (j * DVPAR) * (j * DVPAR);
            }
            dist_DK[i][j] = exp(-(mu_e[i] + 0.5*(j*DVPAR)*(j*DVPAR))) / (2.0 * M_PI);
        }
    }
    for (int i = 0; i < size_mu; i++) {
        for (int j = 0; j < size_vpar; j++) {
            double vpar_mp_j = sqrt(2.0 * U_e_DS[j] - 2.0 * phi_DSE);
            dist[i][j] = bilin_interp(mu_e[i], vpar_mp_j, dist_DK,
                                      mu_e, vpar_MP, size_mu, size_vpar, -1, -1);
        }
    }

    /* output arrays — sized for extended grid */
    int max_op = n_ext + 64;
    double *ne     = calloc(n_ext, sizeof(double));
    double *corr_d = calloc(n_ext, sizeof(double));
    double *corr_c = calloc(n_ext, sizeof(double));
    double *vy_wall = malloc(max_op * sizeof(double));
    double *mu_op   = malloc(max_op * sizeof(double));
    double *chiM    = malloc(max_op * sizeof(double));
    double *dmudvy  = malloc(max_op * sizeof(double));

    mkdir("OUTPUT", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir("TESTS",  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for (int z = 0; z < nzooms; z++) {
        int zoom = zooms[z];
        int size_ne = 0, size_op = 0;
        double flux = 0.0, garbage = 0.0;

        memset(ne,     0, n_ext * sizeof(double));
        memset(corr_d, 0, n_ext * sizeof(double));
        memset(corr_c, 0, n_ext * sizeof(double));

        printf("\n=== zoom = %d ===\n", zoom);
        char mupath[256];
        snprintf(mupath, sizeof(mupath),
                 "OUTPUT/test_DS_mu_%s_zoom%d.txt", ds_base, zoom);
        FILE *fmu = fopen(mupath, "w");
        if (!fmu) fprintf(stderr, "Cannot open %s\n", mupath);
        DENSFINORB(1.0, 1.0, alpha, n_ext, &size_ne,
                   ne, corr_d, corr_c, x_ext, phi_ext, -1.0,
                   dist, mu_e, U_e_DS, size_mu, size_vpar,
                   0.0, &flux, &garbage, zoom,
                   margin, -999.9,
                   vy_wall, mu_op, chiM, dmudvy, &size_op,
                   fmu, NULL);
        if (fmu) { fclose(fmu); printf("Wrote %s\n", mupath); }
        printf("size_ne = %d, flux = %.6f\n", size_ne, flux);

        int n_out = (size_ne < n_phi) ? size_ne : n_phi;
        char outpath[256];
        snprintf(outpath, sizeof(outpath),
                 "OUTPUT/test_DS_ne_%s_zoom%d.txt", ds_base, zoom);
        FILE *fo = fopen(outpath, "w");
        if (!fo) { fprintf(stderr, "Cannot open %s\n", outpath); continue; }
        fprintf(fo, "# x phi ne corr_d corr_c\n");
        for (int i = 0; i < n_out; i++)
            fprintf(fo, "%.12e %.12e %.12e %.12e %.12e\n",
                    x_DS[i], phi_DS[i], ne[i], corr_d[i], corr_c[i]);
        fclose(fo);
        printf("Wrote %s (%d points)\n", outpath, n_out);
    }

    for (int i = 0; i < size_mu; i++) { free(dist[i]); free(dist_DK[i]); }
    free(dist); free(dist_DK); free(mu_e); free(U_e_DS); free(vpar_MP);
    free(x_DS); free(phi_DS);
    free(x_ext); free(phi_ext);
    free(ne); free(corr_d); free(corr_c);
    free(vy_wall); free(mu_op); free(chiM); free(dmudvy);
    free(zooms);
    return 0;
}
