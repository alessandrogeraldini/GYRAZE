# MagSheath

There are currently two different versions of the code: a 1D version called MPS and a 2D version called MPS2d. They are entirely separate, using different scripts and a different executable. This choice was made because the 2D version involves highly nontrivial changes, and so it was decided that it was safer to have an independent 2D version simplified to deal with Boltzmann electrons, and only at a much later stage merge the 2D feutures as an option in a single code. The 2D version is a work-in-progress project, not meant to be used by anyone but the developers. The standard 1D version is written such that it can be run by external users, with the exception of the parameter rho_e/lambda_D (3rd line in input file) which should for the moment be kept to zero by anyone but the developers.

## Input file

The file inputfile.txt contains nine lines:

1st line: a specifier: e.g. ADHOC uses analytical ad hoc distribution functions at the magnetic presheath entrance. If anything other than ADHOC is written here, the code can only run if the files Fi_mpe.txt and Fi_mpe_args.txt are present in the folder. If ADHOC is written here, these files are not necessary because the code defines the ad hoc distribution functions by itself.

2nd line: the angle between the magnetic field angle and the wall, measured in degrees

3rd line: the electron gyroradius rho_e normalised by the Debye length lambda_D

4th line: the number of ion species present

5th line: the ratio of the density of the each species, separated by a space

6th line: the temperature of each species as a fraction of the electron temperature, separated by a space

7th line: the mass of each ion species as a multiple of the electron mass, separated by a space

8th line: ( 1 / 0 ) if ( wall current / electrostatic potential ) is specified on the next line

9th line: the value of current or electrostatic potential

To solve only the magnetic presheath (without the full Debye sheath potential profile) with a simplified electron model it is sufficient to set the value of rho_e/lambda_D to 0.0 or to 100.0 for the two simplified electron models. Setting to 0.0 is more standard.

## Ion distribution function files

These files are only necessary, and opened by the code, if the first line of the input file does not read "ADHOC". For the moment the code can only read a set of distribution function files from one ion species. The change required to read more sets of files corresponding to more ion species are expected to be straightforward.

We have two text files:

Fi_mpe_args.txt includes the grid of values of (1/2) m_i v_perp^2 / T_i in the first line, and the grid of values of (1/2) m_i v_parallel^2 / T_i in the second line. In each line, the values should be separated by a space and there should not be a space after the last value on each line.

Fi_mpe.txt includes the values of the ion distribution function defined on the grids which are read in the Fi_mpe_args.txt file. The number of columns of the file is the number of grid points in (1/2) m_i v_parallel^2 / T_i (second line in Fi_mpe_args), and the number of rows is the number of grid points in (1/2) m_i v_perp^2 / T_i (first line in Fi_mpe_args.txt). Values in each row should be separated only by a space, and there should be no space after the last value and no blank line after the last line of values.

## Electron distribution function files

These files are only necessary, and opened by the code, if the first line of the input file does not read "ADHOC". 

Fe_mpe_args.txt includes the grid of values of (1/2) m_e v_perp^2 / T_e in the first line, and the grid of values of v_parallel / sqrt(T_e/m_e) in the second line. In each line, the values should be separated by a space and there should not be a space after the last value on each line.

Fe_mpe.txt includes the values of the electron distribution function defined on the grids which are read in the Fe_mpe_args.txt file. The number of columns of the file is the number of grid points in v_parallel / sqrt(T_e/m_e) (second line in Fi_mpe_args), and the number of rows is the number of grid points in (1/2) m_e v_perp^2 / T_e (first line in Fe_mpe_args.txt). Values in each row should be separated only by a space, and there should be no space after the last value and no blank line after the last line of values.
