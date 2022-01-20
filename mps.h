#define DEBUG 0
// used to make print statements appear when needed
#define TINY 1e-12
// used to make some inequalities work numerically in case of exact equality.

double tophat(double x1, double x2, double x); 
void newguess(int *convergence, double Te, double *x_grid, double* ne_grid, double *ni, double* phi_grid,int p_size, int *size_ngrid, double lambdaDoverl, double v_cutDS, double pfac, double weight); // gsl_permutation *p, gsl_matrix *m);
void newvcut(double *v_cut, double v_cutDS, double mioverme, double u_i, double u_e, double current, double error_current, int *convergence, double weight);
void denszeroorb(double charge, double *phi,double *n_grid, int p_size, double *Phi_e_point, double **distfunc, double *vpar, double *mu, int size_vpar, int size_mu, double *mue_cut_lookup, double *vpar_cut_lookup, int size_cut);
double *linetodata(char line[], int lenline, int *size);
//double *linetodatanew(char *line, int *size);
void densfinorb(int output_type, int output_dim, double Te, double alpha, int size_phigrid, int *size_ngrid, double* n_grid, double *x_grid, double* phi_grid, double charge, double **FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *vx_i_DS, double *fx_i_DS, double *flux, int zoomfactor, double stopdens);
double bilin_interp(double x, double y, double **FF, double *xx, double *yy, int cols, int rows, int guessi, int guessj);
double lin_interp(double* x_grid, double* y_grid,double given_x ,int n, int line);
