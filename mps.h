#define DEBUG 0
// used to make print statements appear when needed
#define TINY 1e-12
// used to make some inequalities work numerically in case of exact equality.
#define MUGAUSSQUAD 1
// 1: mu of closed orbits in densfinorb(_par) by Gauss-Chebyshev quadrature (mu_gaussquad); 0: original grid integration
#define MUGAUSSN 16
// quadrature points per orbit when MUGAUSSQUAD == 1
#define ANDERSON_M 5
// Picard DS update (ds solver 0): mix the last M iterates (Anderson acceleration); 0 = plain damped Picard
#define NR_SAFESTEP 1
// DS Newton (ds solver 1): 1 = steps accepted/rejected on the true residual, with a trust region; 0 = original frozen-density backtracking
#define DS_DPHIMAX 0.005
// largest change in phi allowed in one DS update (Newton with NR_SAFESTEP 1, or Anderson)
#define NONLOCAL_JAC 12
// DS Newton (ds solver 1): number of localized phi perturbations used to measure the nonlocal
// electron response dne/dphi; 0 = local (tridiagonal) Jacobian, as before
#define NONLOCAL_JAC_EPS 0.001
// amplitude of those perturbations
#define NONLOCAL_JAC_EVERY 1
// reuse a measured response for this many DS Newton steps before measuring it again
#define NONLOCAL_JAC_IONS 1
// with NONLOCAL_JAC > 0: measure the ion response in the same perturbation and build the Jacobian from
// the measured d(ne - ni)/dphi. The two responses nearly cancel, so modelling them separately leaves a
// large relative error in their difference, which is the term the Newton system depends on.
// 0 = measure n_e only and take dni/dphi from ni_corr (analytic)
#define NONLOCAL_JAC_STEPDIR 0
// with NONLOCAL_JAC > 0: also measure the response along the previous Newton direction and correct the
// Jacobian along it (one more density evaluation per step); 0 = bump columns only
#define LM_LAMBDA 0.3
// DS Newton (ds solver 1): Levenberg-Marquardt damping of the linear solve; singular values of the Jacobian
// well above LM_LAMBDA keep the full Newton step, the near-singular smooth mode (0.02-0.16 in the cases
// looked at, next one ~1) is damped. 0 = plain Newton (LU solve)
#define DS_RESTART_CORRECTION 0
// on a restart whose wall potential differs from the new target: 1 = correct the restart phi with a
// linearized BVP solve before the first DS iteration (correct_phi_DS_restart); 0 = start from the restart
// phi unchanged and let the DS solver move the wall gradually
#define DS_WALL_TOL 1e-3
// the DS only counts as converged once its wall potential is within this of -0.5 v_cutDS^2 (<= 0: no check)
#define MUGAUSSSPLINE 0
// phi between grid points in mu_gaussquad: 1 = natural cubic spline, 0 = linear interpolation

struct distfuncDKGK { // contains the distribution function on a 2D grid and the corresponding grid
	double **F;
	double *perp;
	double *par;
	int len_perp;
	int len_par;
}; 

double tophat(double x1, double x2, double x); 
void newguess(double *x_grid, double* ne_grid, double *ni, double* phi_grid,int p_size, int size_ngrid, double lambdaDoverl, double v_cutDS, double pfac, double weight); // gsl_permutation *p, gsl_matrix *m);
void newguess_NR(double *x_grid, double *ne_grid, double *ni_grid, double *phi_grid, int size_phigrid, int size_ngridin, double invgammasq, double v_cutDS, double pfac, double weight, double *ne_corr_delta, double *ne_corr_chiM, double *ni_corr, double *jac_y, double *jac_h, int jac_K, double *dir_h, double *dir_y);
void correct_phi_DS_restart(double *x_DSgrid, double *phi_DSgrid, int size_phiDSgrid, double invgammasq, double v_cutDS, double *ne_DSgrid, double *ne_DSgrid_corr_delta, double *ne_DSgrid_corr_chiM, double *ni_corr, int N_bvp);
void newvcut(double *v_cut, double v_cutDS, double u_i, double u_e, double current, double error_current, double weight);
void error_Poisson(double *error, double *x_grid, double *ne_grid, double *ni_grid, double *nioverne, double *phi_grid, int size_phigrid, int size_ngrid, double invgammasq);
void denszeroorb(double charge, double TeovTs, double *phi,double *n_grid, int p_size, double *Phi_e_point, double *Qe_point, double **distfunc, double *vpar, double *mu, int size_vpar, int size_mu, double *vpar_cut_lookup, double gamma, double *x_grid, double *n_inf); 
double *linetodata(char line[], int lenline, int *size);
//double *linetodatanew(char *line, int *size);
void densfinorb(double Te, double lenfactor, double alpha, int size_phigrid, int *size_ngrid, double* n_grid, double* n_grid_corr_delta, double* n_grid_corr_chiM, double *x_grid, double* phi_grid, double charge, double **FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *flux, double *Qflux, int zoomfactor, double stopdens, double phiDSbump, double *vy, double *mu_op, double *chiMax_op, double *dmudvy, int *size_op, FILE *fmu, FILE *fjmc_out);
void densfinorb_par(double Ti, double lenfactor, double alpha, int size_phigrid, int *size_ngrid, double *n_grid, double *n_grid_corr_delta, double *n_grid_corr_chiM, double *x_grid, double *phi_grid, double charge, double **FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *flux, double *Qflux, int zoomfactor, double margin, double phi_DSbump, double *vy_op, double *mu_op, double *chiMax_op, double *dmudvy_op, int *size_op, FILE *fmu, FILE *fjmc_out);
double bilin_interp(double x, double y, double **FF, double *xx, double *yy, int cols, int rows, int guessi, int guessj);
double lin_interp(double* x_grid, double* y_grid,double given_x ,int n, int line);
void mu_gaussquad(double **mu, double **Uperp, double *xbar, double *xx, double *phi, int size_grid, int sizexbar, int *imax, int *imin, int *kdrop, int *upperlimit, int n);
void densionDS2(double alpha, double TiovTe, double *Bohm, double *ni_DS, double *phi_DS, double phi0, double **FF, double *mu, double *Uminmu, double *vy, double *mu_op, double *chiM, double *twopidmudvy, int size_phi, int size_mu, int size_U, int size_op_i, double* ni_DScorr, double *ni_DS_reflected);
void Figen2(double ***ffarr, double **Uminmuarr, double **muarr, int num_spec, double *nioverne, double *mioverme, double *TioverTe, int *sizevpar, int *sizevperp, double dvpar, double dvperp);
