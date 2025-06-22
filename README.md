# GYRAZE: the grazing-angle gyrokinetic sheath code

This code solves the magnetised SHEATH region, composed of the magnetic presheath, whose thickness is of the order of several ion gyroradii, and the Debye sheath, whose thickness is of the order of several Debye lengths << ion gyroradius. It exploits the asymptotic limit lambda_D / rho_i --> 0 and the GRAZING-ANGLE alpha << 1 of the magnetic field with the wall. 
This means that the Debye sheath and magnetic presheath potential profiles are solved separately, and the ion density in the magnetic presheath and electron density in the Debye sheath are solved using a GYROKINETIC treatment exploiting gyrophase independence and an adiabatic invariant to lowest order in alpha. 
For more details on the theory, see papers by Geraldini, Parra, Militello and by Ewart, Parra and Geraldini. The paper explaining the solution of the combined Debye sheath and magnetic presheath will soon be submitted for publication, and will therefore soon be available on arxiv.com.

GYRAZE can obtain quick solutions for electrostatic potential and density profiles in the magnetised sheath, and can also provide particle distribution functions leaving the magnetised sheath. Ion distribution functions hitting the target are useful for sputtering predictions, while the reflected electron velocity distribution is needed to calculate the electron fluxes through the sheath.

A caveat is that the code cannot capture non-monotonic sheath profiles, which may arise at very shallow angles. However, the code can predict critical angles below which a monotonic solution cannot exist in the magnetised sheath. 
Another caveat is that the code transitions to being very inaccurate at magnetic field angles of 5-8 degrees. 
Above such angles, a solver which does not hard-wire grazing incidence is needed.


## Input file

The file input_physparams.txt contains the physical parameters of the simulation: 1 string and 8 numerical parameters:

string parameter = type_dist_func: this can be "ADHOC" to tell the code to use analytical ad hoc distribution functions at the magnetic presheath entrance. If anything other than ADHOC is written here, the code can only run if the files Fi_mpe.txt and Fi_mpe_args.txt (plus the corresponding electron files starting in Fe) are present in the folder. If ADHOC is written here, these files are not necessary because the code has some built-in ad hoc distribution functions (see paper Geraldini, Parra, Militello 2019).

1st numerical parameter = alphadeg: the angle between the magnetic field angle and the wall, measured in degrees

2nd numerical parameter = gamma_ref: the electron gyroradius rho_e normalised by the Debye length lambda_D

3th numerical parameter = nspec: the number of ion species present

4th numerical parameter = nioverne: the ratio of the density of the each species, separated by a space

5th numerical parameter = TioverTe: the temperature of each species as a fraction of the electron temperature, separated by a space

6th numerical parameter = mioverme: the mass of each ion species as a multiple of the electron mass, separated by a space

7th numerical parameter = set_current: ( 1 / 0 ) if ( wall current / electrostatic potential ) is specified on the next line

8th numerical parameter = target_current/phi_wall: the value of current or electrostatic potential at the wall

To solve only the magnetic presheath (without the full Debye sheath potential profile) with a simplified electron model it is sufficient to set the value of rho_e/lambda_D to 0.0.

The parameter tau must be kept above ~0.2. Simulations crash when tau is too small because the small-alpha equations to lowest order cannot recover the fluid limit (they start missing important terms, see Geraldini, Parra and Militello 2019).

Simulations crash or don't converge if the angle alpha is below the critical angle. For larger values of the parameter gamma the critical angle becomes larger. Therefore, use parameter gamma with caution. Below the critical angle, the code is unable to obtain the solution (because the solution is probably non-monotonic, while the code only handles monotonic profiles).

Other parameter choices that might cause the simulation to crash: a very large electron current, a very large ion current (bounded by the incoming distribution function, which for the moment is independent of the current), a value of wall potential which is not negative enough (the sheath is non-monotonic below a larger critical angle in this case).

Simulations can run with multiple ion species at the moment only with a specific choice of distribution functions. More work is needed to include other options in how the combined distribution functions satisfy the Chodura condition. The output distribution function file is only produced for one ion distribution function.

## Ion distribution function files

These files are only necessary, and opened by the code, if the first line of the input file does not read "ADHOC". For the moment the code can only read a set of distribution function files from one ion species. The change required to read more sets of files corresponding to more ion species are expected to be straightforward.

We have two text files:

Fi_mpe_args.txt includes the grid of values of (1/2) m_i v_perp^2 / T_i in the first line, and the grid of values of (1/2) m_i v_parallel^2 / T_i in the second line. In each line, the values should be separated by a space and there should not be a space after the last value on each line.

Fi_mpe.txt includes the values of the ion distribution function defined on the grids which are read in the Fi_mpe_args.txt file. The number of columns of the file is the number of grid points in (1/2) m_i v_parallel^2 / T_i (second line in Fi_mpe_args), and the number of rows is the number of grid points in (1/2) m_i v_perp^2 / T_i (first line in Fi_mpe_args.txt). Values in each row should be separated only by a space, and there should be no space after the last value and no blank line after the last line of values.

## Electron distribution function files

These files are only necessary, and opened by the code, if the first line of the input file does not read "ADHOC". 

Fe_mpe_args.txt includes the grid of values of (1/2) m_e v_perp^2 / T_e in the first line, and the grid of values of v_parallel / sqrt(T_e/m_e) in the second line. In each line, the values should be separated by a space and there should not be a space after the last value on each line.

Fe_mpe.txt includes the values of the electron distribution function defined on the grids which are read in the Fe_mpe_args.txt file. The number of columns of the file is the number of grid points in v_parallel / sqrt(T_e/m_e) (second line in Fi_mpe_args), and the number of rows is the number of grid points in (1/2) m_e v_perp^2 / T_e (first line in Fe_mpe_args.txt). Values in each row should be separated only by a space, and there should be no space after the last value and no blank line after the last line of values.

## Compiling

Dependencies:
- GNU Scientific Library (GSL): `gsl`, <https://www.gnu.org/software/gsl/>
- CBLAS: `cblas`, <https://www.netlib.org/lapack>

Compile for release (optimization flags)

```sh
make PROFILE=release
```

Compile for debug (requires `libasan`)

```sh
make PROFILE=debug
```
