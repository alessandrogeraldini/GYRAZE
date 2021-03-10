struct spline {
	double* a;
	double* b;
	double* c;
	double* d;

};
double tophat(double x1, double x2, double x); 
int newguess(int convergence, int problem, double* ne_grid, double* phi_grid,int p_size,double* ne_new_grid, double * phi_new_grid, double* ne_cut_grid, double* phi_cut_grid, int cut_points,int *use_new, double u_e, double *v_cut);
void makelookup(double *phi,double *ne, int p_size,double v_cut, double v_cutDS, double *Phi_e_point);
double *linetodata(char line[], int lenline, int *size);
int densfinorb(int dodistfunc,double* ne_grid, double* phi_grid, int *p_size);
void cutlookup(double* phi_cut, double* ne_cut, int cut_points);
void Fgenerator(double phi_cut);
double bilin_interp(double x, double y, double **FF, double *xx, double *yy, int cols, int rows, int guessi, int guessj);
void generate_cspline(struct spline *spl,double *x,double *y, int n);
double lin_interp(double* x_grid, double* y_grid,double given_x ,int n, int line);
