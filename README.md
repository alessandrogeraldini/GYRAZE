# MagSheath

There are currently 3 different versions of the code: an old (discontinued) version called MPS which simply had different normalizations in the input and output of the scripts evaluating the ion density in the magnetic presheath and the electron density in the Debye sheath (will soon be deleted); a newer version called MPS_renorm; a 2d (in space) version called MPS2d.

## Input file
The file inputfile.txt contains seven lines: angle alpha in degrees, electron gyroradius over Debye length, number of ion species, ion to electron temperature, ion to electron mass, a constant equal to 1 when setting current to wall or 0 wheh setting the wall potential, the value of current or potential (whichever was indicated in the previous line).

Lines in the input file:

1st line: a specifier: e.g. ADHOC uses analytical ad hoc distribution functions at the magnetic presheath entrance. If anything other than ADHOC is written here, the code can only run if the files Fi_mpe.txt and Fi_mpe_args.txt are present in the folder. If ADHOC is written here, these files are not necessary because the code defines the adhoc distribution functions by itself.

2nd line: the angle between the magnetic field angle and the wall, measured in degrees

3rd line: the electron gyroradius rho\_e normalised by the Debye length lambda_\D

4th line: the number of ion species present

5th line: the ratio of the density of the each species, separated by a space

6th line: the temperature of each species as a fraction of the electron temperature, separated by a space

7th line: the mass of each ion species as a multiple of the electron mass, separated by a space

8th line: ( 1 / 0 ) if ( wall current / electrostatic potential ) is specified on the next line

9th line: the value of current or electrostatic potential

To solve only the magnetic presheath (without the full Debye sheath potential profile) with a simplified electron model it is sufficient to set the value of gamma = 0.0 or gamma = 100.0 (or larger)
