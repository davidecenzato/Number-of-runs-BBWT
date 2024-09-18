# BBWT_paper

This is the code for the paper **On the Number of Equal-Letter Runs of the Bijective Burrows-Wheeler
Transform** by Elena Biagi, Davide Cenzato, Zsuzsanna Lipták and Giuseppe Romana.

(TODO) brief intro

This code investigates:
- $\rho$: (TODO) runs ratio, formula here
* $\delta$: (TODO) difference in the number of runs of the BBWT formula here
+ $\matchcal{L}$ /# Lyndon factors (TODO fix) $\matchcal{l}(s)-\matchcal{l}(s^{rev})$

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
- Statistics on the runs ratio for each string length up to mak_k (max $\rho$, avg $\rho$, std $\rho$, % strings with $\rho$=1)
* Strings with max $\rho$ for each length
+ Strings with max $\rho$ over all the string up to length max_k
- Statistics on $\delta$ for each string length up to mak_k (max $\delta$, min $\delta$, avg $\delta$, std $\delta$, % strings with $\delta$=0)
* Strings with max $\delta$ for each length
+ Strings with max $\delta$ over all the string up to length max_k
- Strings with min $\delta$ for each length
* Strings with min $\delta$ over all the string up to length max_k

- Statistics on $\matchcal{L}$ for each string length up to mak_k (max $\matchcal{L}$, min $\matchcal{L}$, avg $\matchcal{L}$, std $\matchcal{L}$)
* Strings with max $\matchcal{L}$ for each length
+ Strings with max $\matchcal{L}$ over all the string up to length max_k
- Strings with min $\matchcal{L}$ for each length
* Strings with min $\matchcal{L}$ over all the string up to length max_k
* (TODO add distinc Lyndon factors)
+ (TODO make thic coherent with real results)

Example:
``` console
python BBWT_tests.py 2 6 -o outdir/
```
This will create a subdirectory inside outdir ``` outdir/sigma_2_6 ``` containing the output files for binary strings of lenght from 3 to 6. 