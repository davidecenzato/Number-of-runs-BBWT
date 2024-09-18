# BBWT_paper

This is the code for the paper **On the Number of Equal-Letter Runs of the Bijective Burrows-Wheeler
Transform** by Elena Biagi, Davide Cenzato, Zsuzsanna Lipták and Giuseppe Romana.

(TODO) brief intro

This code investigates:
- $\rho_B$: (TODO) runs ratio, formula here
* $\delta_B$: (TODO) difference in the number of runs of the BBWT formula here
+ $\mathcal{L}$ /# Lyndon factors (TODO fix) $\mathcal{l}(s)-\mathcal{l}(s^{rev})$

### Install repo

```console
git clone https://github.com/davidecenzato/BBWT_paper.git
git submodule update --init --recursive

cd cais

cd external/sdsl-lite
mkdir installed
./install.sh installed/

cd ../../
cp ../main.cpp ./
cp ../Makefile ./
make 

cd ..
mkdir outdir
```

### Running experiments
The code takes in input the size of the alphabet ($\sigma$) and the maximumg length of the strings we are interested in (max_k).
It also creates a temporaty directory that will then be deleted automatically.
```
usage: BBWT_tests.py [-h] [-o OUTDIR] sigma max_k

positional arguments:
  sigma                 size of the alphabet
  max_k                 maximum strings length

options:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        existing output directory for generated files
```

The outputs will be found in a new directory. It consists of:
- Statistics on the runs ratio for each string length up to mak_k (max $\rho_B$, avg $\rho_B$, std $\rho_B$, % strings with $\rho_B$=1)
* Strings with max $\rho_B$ for each length
+ Strings with max $\rho_B$ over all the string up to length max_k
- Statistics on $\delta_B$ for each string length up to mak_k (max $\delta_B$, min $\delta_B$, avg $\delta_B$, std $\delta_B$, % strings with $\delta_B$=0)
* Strings with max $\delta_B$ for each length
+ Strings with max $\delta_B$ over all the string up to length max_k
- Strings with min $\delta_B$ for each length
* Strings with min $\delta_B$ over all the string up to length max_k

- Statistics on $\mathcal{L}$ for each string length up to mak_k (max $\mathcal{L}$, min $\mathcal{L}$, avg $\mathcal{L}$, std $\mathcal{L}$)
* Strings with max $\mathcal{L}$ for each length
+ Strings with max $\mathcal{L}$ over all the string up to length max_k
- Strings with min $\mathcal{L}$ for each length
* Strings with min $\mathcal{L}$ over all the string up to length max_k
* (TODO add distinc Lyndon factors)
+ (TODO make thic coherent with real results)

Example:
``` console
python BBWT_tests.py 2 6 -o outdir/
```
This will create a subdirectory inside outdir ``` outdir/sigma_2_6 ``` containing the output files for binary strings of lenght from 3 to 6. 