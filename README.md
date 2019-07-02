# Files in the ratemeasure repository
The file ratemeasure.py contains Python code to calculate a rate constant estimate based on a least squares optimization for a given mechanism and a given set of experimental conditions. The Mathematica notebooks combustion.nb and pyrolysis.nb contain code for plotting results. The folder experiments contains files with experimental conditions, and the folder mechanisms contains Cantera mechanisms.

# System requirements
The python code has been run with anaconda, which can be downloaded here: https://www.anaconda.com/distribution/. The script requires packes numpy, scipy and cantera, which can be installed in a new environment after installing anaconda with the shell command  
`conda create -n cantera_env -c cantera numpy scipy cantera`  
Next, activate the environment with  
`source activate cantera_env`  

# Usage
Running `./ratemeasure.py -h` produces the following usage message:
```
usage: ratemeasure.py [-h] --filebase FILEBASE [--mechanism MECHANISM]
                      [--Npoints NPOINTS] [--experiments EXPERIMENTS]
                      [--measure MEASURE] [--sensitivity SENSITIVITY]
                      [--remove REMOVE] [--retain RETAIN] [--seed SEED]
                      [--ktol KTOL] [--out {0,1}] [--maxes MAXES [MAXES ...]]
                      [--yields YIELDS [YIELDS ...]]

Measure a rate constant in an incomplete mechanism.

optional arguments:
  -h, --help            show this help message and exit
  --filebase FILEBASE   String for the output files.
  --mechanism MECHANISM
                        Cantera mechanism file. Default mechanisms/gri30.cti.
  --Npoints NPOINTS     Number of time points to output. Default 5e3.
  --experiments EXPERIMENTS
                        File containing a line [tmax temperature pressure
                        initials] for each experimental condition. Initials is
                        a list of pairs [index] [proportion] specifying
                        species with nonzero initial mole fraction. Default
                        experiments/air.dat.
  --measure MEASURE     the index of the reaction whose rate is measured.
                        Default 37.
  --sensitivity SENSITIVITY
                        the index of the reaction whose rate is measured.
                        Default 2.
  --remove REMOVE       The number of reactions to randomly remove. Default
                        40.
  --retain RETAIN       The number of reactions to retain. Default 40.
  --seed SEED           The random seed. Default 1.
  --ktol KTOL           The tolerance in k relative to k0 in the minimization.
                        Default 1e-6.
  --out {0,1}           Flag for outputting: 1 for observation and fit output.
                        Default 1.
  --maxes MAXES [MAXES ...]
                        Indices of concentration maxima to include in error.
                        -1 for none. Default 2.
  --yields YIELDS [YIELDS ...]
                        Indices of concentration yields to include in error.
                        -1 for none. Default -1.
  ```
  
After running once with sensitivities, the files filebasemi.dat and filebasems.dat are created to store sensitivity information.  If these files are present, the program runs without calculating sensitivities to speed up evaluation.  The times, species names, molar fractions, reaction rates, products, reactants, temperatures, and pressures are stores in files based on filebase.

  -----------
# Examples
To find the rate estimate for reaction 37 of the GRI mechanism removing 40 reactions, retaining 40 reactions, using [O] sensitivity, with [O] maxima and no yields, run

`./ratemeasure.py --mechanism mechanisms/gri30.cti --experiments experiments/air.dat --yields -1 --maxes 2 --sensitivity 2 --filebase data/air/ --retain 40 --remove 40 --measure 37`  

To find the rate estimate for reaction 478 of the pyrolysis mechanism removing 200 reations, retaining 200 reactions, using H2 sensitivity, with no maxima and yields of H2 and CH4, run

`./ratemeasure.py --mechanism mechanisms/gri30.cti --experiments experiments/pyrolysis.dat --yields 4 10 --maxes -1 --sensitivity 4 --filebase data/pyrolysis/ --retain 200 --remove 200 --measure 478`  

