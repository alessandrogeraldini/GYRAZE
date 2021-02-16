struct spline {
	double* a;
	double* b;
	double* c;
	double* d;

};
void generate_cspline(struct spline *spl,double *x,double *y, int n);
double lin_interp(double* x_grid, double* y_grid,double given_x ,int n, int line);
