# Study on runs of the BBWT

This is the code for the paper **On the Number of Equal-Letter Runs of the Bijective Burrows-Wheeler
Transform** by Elena Biagi, Davide Cenzato, Zsuzsanna Lipták and Giuseppe Romana.

(TODO) brief intro, ex $r_B$

This code investigates:
- BBWT runs ratio: $\rho_B(s) = \max(\frac{r_B(s)}{r_B(s^{rev})}, \frac{r_B(s^{rev})}{r_B(s)})$
* BBWT runs difference: $\delta_B(s) = r_B(s) - r_B(s^{rev})$
+ Lyndon factors difference: $\ell(s)-\ell(s^{rev})$
- Distinct Lyndon factors difference: $d(s)-d(s^{rev})$

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
- Statistics on the runs ratio for each string length up to mak_k (max, min, avg, std, % strings with $\rho_B(s)=1$)
* Strings with max $\rho_B(s)$ for each length
+ Strings with max $\rho_B(s)$ over all the string up to length max_k
- Statistics on $\delta_B(s)$ for each string length up to mak_k (max, min, avg, std, % strings with $\delta_B(s)=0$)
* Strings with max $\delta_B(s)$ for each length
+ Strings with max $\delta_B(s)$ over all the string up to length max_k
- Strings with min $\delta_B(s)$ for each length
* Strings with min $\delta_B(s)$ over all the string up to length max_k

- Statistics on $\ell(s)-\ell(s^{rev})$ for each string length up to mak_k (max, min, avg, std, % strings with $\ell(s)-\ell(s^{rev})=0$)
* Strings with max $\ell(s)-\ell(s^{rev})$ for each length
+ Strings with max $\ell(s)-\ell(s^{rev})$ over all the string up to length max_k
- Strings with min $\ell(s)-\ell(s^{rev})$ for each length
* Strings with min $\ell(s)-\ell(s^{rev})$ over all the string up to length max_k

- Statistics on $d(s)-d(^{rev})$ for each string length up to mak_k (max, min, avg, std, % strings with $d(s)-d(^{rev})=0$)
* Strings with max $d(s)-d(s^{rev})$ for each length
+ Strings with max $d(s)-d(s^{rev})$ over all the string up to length max_k
- Strings with min $d(s)-d(s^{rev})$ for each length
* Strings with min $d(s)-d(s^{rev})$ over all the string up to length max_k

+ (TODO make this coherent with real results)

Example:
``` console
python BBWT_tests.py 2 6 -o outdir/
```
This will create a subdirectory inside outdir ``` outdir/sigma_2_6 ``` containing the output files for binary strings of lenght from 3 to 6. 