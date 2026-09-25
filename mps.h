#define DEBUG 0
// used to make print statements appear when needed
#define TINY 1e-12
// used to make some inequalities work numerically in case of exact equality.
#ifndef DVZ_QUAD
#define DVZ_QUAD 0.1
#endif
// step in v_z of the trapezoid integrals over the distribution in densfinorb(_par) (closed and open orbits, n_inf,
// electrons and MP ions). Separate from DVPAR, which only sets the grid F is tabulated on (and interpolated from).
// Was hard-coded 0.1. Overridable with -DDVZ_QUAD=... for resolution tests.
#define MUGAUSSQUAD 1
// 1: mu of closed orbits in densfinorb(_par) by Gauss-Chebyshev quadrature (mu_gaussquad); 0: original grid integration
#ifndef MUGAUSSN
#define MUGAUSSN 16
#endif
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
#define NONLOCAL_JAC_BASIS 0
// shape of the NONLOCAL_JAC perturbations. 0 = raised-cosine bumps on disjoint blocks (localized).
// 1 = smooth sine modes sin(pi*(k+1)*j/N) spanning the whole n_e grid, vanishing at the wall and at the
// grid end. Both bases are orthogonal, which is what the rank-1 corrections below need. The disjoint
// bumps narrow as NONLOCAL_JAC grows, so the response per column shrinks while the density noise does
// not, and the Jacobian's smallest singular value sits at that noise floor (~0.06 measured) -- a mode
// indistinguishable from measurement error, which a probe showed the solver cannot actually move along.
// The sine modes measure that smooth subspace directly, at full amplitude, for the same cost per column.
#define NONLOCAL_JAC_EPS 0.002
#define NONLOCAL_JAC_ANALYTIC 1
// DS Newton: build d(ne - ni)/dphi in full instead of the local model plus NONLOCAL_JAC measured directions.
// n_e comes from densfinorb_par itself (dfo_jac): the local chi(x) term, the orbits' mu (an integral of phi
// over the orbit), the accessibility lift of the energy bound (set by a neighbouring orbit's mu, the term
// behind the anti-screening seen at the bump) and the reflection threshold Ucrit(mu), at ~1.5x the cost of
// one density call. n_i by central differences of densionDS (NONLOCAL_JAC_EPS), columns in parallel.
// Checked with test_jac_DS against finite differences of densfinorb_par (3.2, iteration 23): J h within
// 4-10% for sine modes and 5-17% for bumps, where the local model was off by 10-43% and 10-31%. Not
// included: the open-orbit density (x < ~1, where n_e also jumps at x ~ 0.37). Needs use_parallel 1.
// Replaces NONLOCAL_JAC (not measured when this is 1). 0 = off.
// amplitude of those perturbations
#define NONLOCAL_JAC_EVERY 1
// reuse a measured response for this many DS Newton steps before measuring it again
#define NONLOCAL_JAC_IONS 1
// with NONLOCAL_JAC > 0: measure the ion response in the same perturbation and build the Jacobian from
// the measured d(ne - ni)/dphi. The two responses nearly cancel, so modelling them separately leaves a
// large relative error in their difference, which is the term the Newton system depends on.
// 0 = measure n_e only and take dni/dphi from ni_corr (analytic)
#define NONLOCAL_JAC_CHECK 0.2
// with NONLOCAL_JAC > 0: also measure each column at -eps and average the two sides (central difference);
// where the electron +eps and -eps responses differ by more than this, n_e is not smooth in phi there and
// that entry falls back to the local model. 0 = one-sided measurement, no check. The electron response has a
// noise floor ~0.03 at eps = 0.001 (median |y+ - y-|, 3.3 case); 0.2 catches only the real jumps (x = 0.37 near
// the wall, and the n_e grid edge)
#define NONLOCAL_JAC_STEPDIR 0
// with NONLOCAL_JAC > 0: also measure the response along the previous Newton direction and correct the
// Jacobian along it (one more density evaluation per step); 0 = bump columns only
#define LM_LAMBDA 0.3
// DS Newton (ds solver 1): Levenberg-Marquardt damping of the linear solve; singular values of the Jacobian
// well above LM_LAMBDA keep the full Newton step, the near-singular smooth mode (0.02-0.16 in the cases
// looked at, next one ~1) is damped. 0 = plain Newton (LU solve)
#define LM_LAMBDA_END 0.0
// LM damping used once the wall is within DS_WALL_TOL of its target (the continuation is then over and what
// is left lives in the near-singular mode): smaller than LM_LAMBDA lets that mode move. 0 = keep LM_LAMBDA
#define DS_LINESEARCH 0
// 1 = after each DS Newton step, search along the Jacobian's weakest singular vector using true residual
// 2 = diagnostic: probe that mode out to the amplitude the linear model asks for and print the whole
// curve. The model claims a step of |u^T F|/s there removes most of the residual while a probe at 0.01
// made it worse; only an evaluation at the model's own amplitude can say which is right. Original note:
// after each DS Newton step, search along the Jacobian's weakest singular vector using true residual
// evaluations. LM damps that mode by ~s/lambda^2 relative to Newton, and at a stall s falls to ~1e-3, so
// the step cannot traverse it; the amplitude is set by evaluating the residual rather than by the
// linearization, whose measured d(ne-ni)/dphi is least reliable in exactly that mode. 0 = off
#define DS_LS_T0 0.01
// first amplitude probed, as max |dphi| along the mode (deliberately above DS_DPHIMAX: the point is to
// move further along this one direction than the trust region allows, with an evaluation to justify it)
#define DS_LS_TMAX 0.08
// largest amplitude the search may reach by doubling
#define DS_LS_TAPER 0
// number of grid points over which the weak-mode vector is tapered to zero before the n_e grid end, so a
// step along it leaves the asymptotic tail alone. 0 = no taper. The taper injects curvature of order
// amplitude/(DS_LS_TAPER*dx)^2 right at the edge, which is a residual the probe would then measure
// instead of the mode's real effect -- so the diagnostic probe runs untapered.
#define DS_LS_EDGE_SKIP 3
// grid points excluded on the inner side of the join, where the weak-mode vector falls to zero and the
// perturbed phi meets the unchanged asymptotic tail. That kink dominates the max however small the step
// (|t| = 5e-4 put it at x = 5.93, 0.0228, against 0.0141 at x = 4.17), so the probe measures the mode's
// effect below it. Counted from the last nonzero entry of the mode, NOT from the end of the n_e grid:
// the grid is resized by DENSFINORB under perturbation and the join does not sit at a fixed offset.
#define DS_LS_PROBE_CAP 0.05
// largest |t| the diagnostic probe (DS_LINESEARCH 2) will try. The orbit integration aborts with
// "more than one maximum" somewhere between 0.098 and 0.196, so the model amplitude itself (~0.1) is
// not always reachable; the probe is capped and reports how far it got.
#define DS_LS_SKIP 3
// diagnostic probe: skip this many DS Newton calls first, so it measures the plateau rather than the
// restart transient. The earlier probe measured the transient by accident and its weak mode is a
// different, noise-dominated one (s ~ 0.06) from the plateau's (s ~ 0.17, unchanged by NONLOCAL_JAC_EPS)
#define DS_LS_TMIN 1e-4
// smallest amplitude probed: if neither sign helps at DS_LS_T0 the scale is quartered down to this before
// giving up, since the useful amplitude may be well below DS_DPHIMAX rather than above it
#define DS_LS_MAXEVAL 9
// most density evaluations spent on the search per Newton step (1 baseline + probes)
#define DS_RESTART_CORRECTION 0
// on a restart whose wall potential differs from the new target: 1 = correct the restart phi with a
// linearized BVP solve before the first DS iteration (correct_phi_DS_restart); 0 = start from the restart
// phi unchanged and let the DS solver move the wall gradually
#define DS_MINIMAX 0
// DS Newton: target the max residual, which is what convergence is tested on, instead of the rms.
// (a) Before the damped solve, weight row i by sqrt(w_i), w_i = max(DS_MINIMAX_WFLOOR, (|r_i|/max|r|)^(P-2)),
//     so the LM damping suppresses directions that only serve rows far from the max rather than the ones
//     fixing it. The system is square, so this changes the step only through the damping (an undamped
//     Newton step is weighting-invariant).
// (b) Accept/reject on max|r| instead of the rms, with DS_MINIMAX_MERIT_TOL.
// Solver-free runs at 3.25 (MP fixed) found the plateau is a least-squares optimum (a Gauss-Newton step
// moved the rms 0.2%) while a max-targeted fit lowered max|r| 0.0139 -> 0.0132 at the cost of +5% rms,
// which the rms merit test (DS_MERIT_TOL 3%) would reject. 0 = rms objective as before.
#define DS_SLP 0
// DS Newton endgame (wall within DS_WALL_TOL of its target): take the step from a linear program that
// minimises the max of the linearised relative Poisson residual, r_i(d) = (-F_i + (J d)_i) / (gamma^2 ni),
// over a box trust region, with the mean over error_Poisson's rows capped and a cap on second differences
// of the step. This targets the convergence test (the max) directly; LM minimises a 2-norm, and in this
// square system row weights change the step only through the damping (DS_MINIMAX, which did not work).
// Solver-free, the same LP took max|r| 0.0139 -> 0.0115 in one step at 3.25. Needs GLPK (-lglpk).
// The wall approach still uses LM. Merit test on the max (DS_MINIMAX_MERIT_TOL). 0 = off.
#define DS_SLP_DELTA 0.002
// trust region, max |dphi| per step (~0.0023 is where the linearisation was seen to fail); scaled by the
// accept/reject logic's alpha/weight, so a rejected step halves it
#define DS_SLP_TAPER 1.0
// the bound tapers linearly to 0 over this distance before the last unknown, so the tail join stays put
#define DS_SLP_SMOOTH 2e-4
// cap on |d_{j-1} - 2 d_j + d_{j+1}|: forbids grid-scale zigzag (the needed change is ~2e-5)
#define DS_SLP_ENGAGE 0.04
// the LP takes over only once max|r| is below this (and the wall within DS_WALL_TOL of its target); before
// that the steps are LM. On a cold start at 3.2 it engaged at max|r| ~ 2 and built grid-scale zigzag in phi''.
#define DS_SLP_ZIGZAG 1e-4
// the LP step may not grow the fourth difference of phi (phi'' node-to-node zigzag) at any point beyond
// max(its current size, this); 1e-4 is ~0.012 in phi'' at dx = 0.093, above a smooth solution's ~7e-5
#define DS_SLP_TMAX 0.0115
// the LP first finds the smallest achievable max t*, then minimises the mean with the max held at
// max(t*, DS_SLP_TMAX): set just below tol_DS[1] (0.0125) to leave room for the linearization error. The
// fraction of the pending wall move is an LP unknown too, kept so the remaining gap is below DS_WALL_TOL/2.
#if DS_SLP && !NR_SAFESTEP
#error "DS_SLP needs NR_SAFESTEP 1 (it uses the accept/reject logic for its trust region)"
#endif
#if NONLOCAL_JAC_ANALYTIC && DS_LINESEARCH
#error "DS_LINESEARCH uses the NONLOCAL_JAC measurement set-up, which NONLOCAL_JAC_ANALYTIC skips"
#endif
#define DS_MINIMAX_P 10.0
#define DS_MINIMAX_WFLOOR 1e-3
#define DS_MINIMAX_MERIT_TOL 0.005
// fractional rise in max|r| rejected. The max sits in the interior (x ~ 4), so the n_e-grid resize jitter
// that DS_MERIT_TOL allows for in the rms does not reach it; run-to-run scatter is ~1e-5.
#define DS_MERIT_TOL 0.03
// fractional rise in the rms residual that counts as a step worth rejecting. Not zero: the rms is taken
// over the n_e grid, whose size moves by a point or so each iteration, and changing the sample shifts it
// by a few tenths of a percent on its own. Measured scale: run-to-run scatter ~1e-5, resize jitter ~0.7%,
// a genuinely bad step (lambda = 0.1 at the 3.25 plateau) +176%. 3% sits well clear of both ends.
#ifndef DS_XEND
#define DS_XEND 5.2
#endif
// End the DS n_e grid at the first point with x >= DS_XEND, instead of where n_e first comes within MARGIN_DS
// of n_inf on the rising branch. 0 = the density threshold, as before. That threshold sits on a knife edge:
// near its natural end the deficit 1 - n_e/n_inf is ~0.0019 (tail phi ~ -1.3e-3 plus loss-cone depletion), and
// n_inf can itself sit ~0.2% above the far-field density, so a 1e-4 shift at a restart moved the end from
// x ~ 6.5 to x ~ 16.6, and in steady state the grid length locked into the bump/dip limit cycle (67 <-> 71
// points at 3.25). x = 9 is about where MARGIN_DS = 0.0005 stops on the 3.25 plateau. Handed to densfinorb(_par)
// through its margin argument: a margin >= 1 means "stop at this x" (a negative one still means "don't stop").
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
void newguess_NR(double *x_grid, double *ne_grid, double *ni_grid, double *phi_grid, int size_phigrid, int size_ngridin, double invgammasq, double v_cutDS, double pfac, double weight, double *ne_corr_delta, double *ne_corr_chiM, double *ni_corr, double *jac_y, double *jac_h, int jac_K, double *dir_h, double *dir_y, double *ls_v, double *ls_tmodel, double *jac_e, double *jac_i, int jac_n);
void correct_phi_DS_restart(double *x_DSgrid, double *phi_DSgrid, int size_phiDSgrid, double invgammasq, double v_cutDS, double *ne_DSgrid, double *ne_DSgrid_corr_delta, double *ne_DSgrid_corr_chiM, double *ni_corr, int N_bvp);
void newvcut(double *v_cut, double v_cutDS, double u_i, double u_e, double current, double error_current, double weight);
extern int error_Poisson_imax;   /* grid index of the largest residual from the last error_Poisson call */
void error_Poisson(double *error, double *x_grid, double *ne_grid, double *ni_grid, double *nioverne, double *phi_grid, int size_phigrid, int size_ngrid, double invgammasq);
void denszeroorb(double charge, double TeovTs, double *phi,double *n_grid, int p_size, double *Phi_e_point, double *Qe_point, double **distfunc, double *vpar, double *mu, int size_vpar, int size_mu, double *vpar_cut_lookup, double gamma, double *x_grid, double *n_inf); 
double *linetodata(char line[], int lenline, int *size);
//double *linetodatanew(char *line, int *size);
void densfinorb(double Te, double lenfactor, double alpha, int size_phigrid, int *size_ngrid, double* n_grid, double* n_grid_corr_delta, double* n_grid_corr_chiM, double *x_grid, double* phi_grid, double charge, double **FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *flux, double *Qflux, int zoomfactor, double stopdens, double phiDSbump, double *vy, double *mu_op, double *chiMax_op, double *dmudvy, int *size_op, FILE *fmu, FILE *fjmc_out);
void densfinorb_par(double Ti, double lenfactor, double alpha, int size_phigrid, int *size_ngrid, double *n_grid, double *n_grid_corr_delta, double *n_grid_corr_chiM, double *x_grid, double *phi_grid, double charge, double **FF, double *mumu, double *UU, int sizemumu, int sizeUU, double grid_parameter, double *flux, double *Qflux, int zoomfactor, double margin, double phi_DSbump, double *vy_op, double *mu_op, double *chiMax_op, double *dmudvy_op, int *size_op, FILE *fmu, FILE *fjmc_out);
double bilin_interp(double x, double y, double **FF, double *xx, double *yy, int cols, int rows, int guessi, int guessj);
/* Analytic electron Jacobian from densfinorb_par (see NONLOCAL_JAC_ANALYTIC). When dfo_jac != NULL, the next
 * electron (charge < 0, zoomfactor 1) call also fills dfo_jac[i*dfo_jac_n + m] = d n(x_i)/d phi_grid[m], normalized
 * like n_grid, for the rows it computes (i < *size_ngrid) and columns m < dfo_jac_n. dfo_jac_terms selects the
 * channels (DFO_JAC_* bits; for testing). */
extern double *dfo_jac;
extern int dfo_jac_n, dfo_jac_terms;
#define DFO_JAC_LOCAL   1   /* chi(x) at the point itself: the v_x <-> U_perp map, incl. the chi_M limit there */
#define DFO_JAC_CHIM    2   /* chi_M of each orbit: the upper U_perp limit, delta psi at the orbit's maximum */
#define DFO_JAC_MU      4   /* mu of each closed orbit, an integral of psi over the orbit */
#define DFO_JAC_BARRIER 8   /* the reflection threshold U_crit(mu), set by the open orbits' chi_M and mu */
#define DFO_JAC_OPEN   16   /* the open-orbit density */
#define DFO_JAC_ALL    31
double lin_interp(double* x_grid, double* y_grid,double given_x ,int n, int line);
void mu_gaussquad(double **mu, double **Uperp, double *xbar, double *xx, double *phi, int size_grid, int sizexbar, int *imax, int *imin, int *kdrop, int *upperlimit, int n);
void densionDS2(double alpha, double TiovTe, double *Bohm, double *ni_DS, double *phi_DS, double phi0, double **FF, double *mu, double *Uminmu, double *vy, double *mu_op, double *chiM, double *twopidmudvy, int size_phi, int size_mu, int size_U, int size_op_i, double* ni_DScorr, double *ni_DS_reflected);
void Figen2(double ***ffarr, double **Uminmuarr, double **muarr, int num_spec, double *nioverne, double *mioverme, double *TioverTe, int *sizevpar, int *sizevperp, double dvpar, double dvperp);
