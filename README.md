# MagSheath
## Before running
It is important to first generate the electron and ion distribution functions to be used by the code.
The scripts Fi_mpe_generate.py and Fe_mpe_generate.py take care of this, reading the temperature in the input file.

## Input file
The file inputfile.txt contains six lines: angle alpha in degrees, electron gyroradius over Debye length, number of ion species, ion to electron temperature, ion to electron mass, a constant equal to 1 when setting current to wall or 0 wheh setting the wall potential, the value of current or potential (whichever was indicated in the previous line).

Example input file:

4.0 (angle in degrees)

0.8 (electron gyroradius to Debye length)

1 (number of species, only 1 supported for now)

1.0 (ion to electron temperature ratio)

3600.0 (ion to electron mass ratio)

1 (set current)

0.0 (value of current, in this case zero current = ambipolarity)

To solve only the magnetic presheath (without the full Debye sheath potential profile) with a simplified electron model it is sufficient to set the value of gamma = 0.0 or gamma = 10.0 (or larger)
