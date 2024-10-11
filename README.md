# Study on the number of runs of the BBWT

This is the code to replicate the experiments in the paper **On the Number of Equal-Letter Runs of the Bijective Burrows-Wheeler
Transform** by Elena Biagi, Davide Cenzato, Zsuzsanna Lipták and Giuseppe Romana.

Here, we study the relationship between the number of runs of the BBWT of a text $r_B(s)$ and its reverse $r_B(s^{rev})$ by computing different parameters, including:
- BBWT runs ratio: $\rho_B(s) = \max(\frac{r_B(s)}{r_B(s^{rev})}, \frac{r_B(s^{rev})}{r_B(s)})$
* BBWT runs difference: $\delta_B(s) = r_B(s) - r_B(s^{rev})$
+ Lyndon factors difference: $\ell(s)-\ell(s^{rev})$
- Distinct Lyndon factors difference: $d(s)-d(s^{rev})$

### Download and Compile

```console
git clone https://github.com/davidecenzato/Number-of-runs-BBWT.git
git submodule update --init --recursive

python compile.py
```

### Requirements

This software requires:
* A modern 64-bit Mac OS or Linux operating system.
* A modern Python 3 release version 3.7 or higher.
* A modern C++11 compiler such as `g++` version 4.9 or higher.

### Running experiments
The software takes in input the size of the alphabet ($\sigma$) and the maximum length of the strings we are interested in ($max_k$).
It also creates a temporary directory that will then be deleted automatically.
```
usage: BBWT_tests.py [-h] [-o OUTDIR] sigma max_k

positional arguments:
  sigma                 Size of the alphabet
  max_k                 Maximum string length

options:
  -h, --help            Show this help message and exit
  -o OUTDIR, --outdir OUTDIR
                        Existing output directory for generated files
```
Example:
``` console
python BBWT_tests.py 2 6 -o outdir/
```
This command will create a subdirectory inside /outdir named /outdir/sigma_2_6 containing the output files for binary strings of length 3 to 6.


The output files consist of:
- Statistics on the runs ratio for each string length up to mak_k (max, min, avg, std, % of strings with $\rho_B(s)=1$)
* Strings with max $\rho_B(s)$ for each length
+ Strings with max $\rho_B(s)$ over all the strings up to length $max_k$
- Statistics on $\delta_B(s)$ for each string length up to mak_k (max, min, avg, std, % of strings with $\delta_B(s)=0$)
* Strings with max $\delta_B(s)$ for each length
+ Strings with max $\delta_B(s)$ over all the strings up to length $max_k$
- Strings with min $\delta_B(s)$ for each length
* Strings with min $\delta_B(s)$ over all the strings up to length $max_k$

- Statistics on $\ell(s)-\ell(s^{rev})$ for each string length up to mak_k (max, min, avg, std, % of strings with $\ell(s)-\ell(s^{rev})=0$)
* Strings with max $\ell(s)-\ell(s^{rev})$ for each length
+ Strings with max $\ell(s)-\ell(s^{rev})$ over all the strings up to length $max_k$
- Strings with min $\ell(s)-\ell(s^{rev})$ for each length
* Strings with min $\ell(s)-\ell(s^{rev})$ over all the strings up to length $max_k$

- Statistics on $d(s)-d(^{rev})$ for each string length up to mak_k (max, min, avg, std, % of strings with $d(s)-d(s^{rev})=0$)
* Strings with max $d(s)-d(s^{rev})$ for each length
+ Strings with max $d(s)-d(s^{rev})$ over all the strings up to length $max_k$
- Strings with min $d(s)-d(s^{rev})$ for each length
* Strings with min $d(s)-d(s^{rev})$ over all the strings up to length $max_k$

### References and citations

[1] Elena Biagi, Davide Cenzato, Zsuzsanna Lipták, Giuseppe Romana: On the Number of Equal-Letter Runs of the Bijective Burrows-Wheeler Transform. ICTCS 2023: 129-142 ([go to the paper](https://ceur-ws.org/Vol-3587/4564.pdf))

If you use this work, please cite the following papers:

#### conference version
    @inproceedings{BiagiCLR23,
      author       = {Elena Biagi and
                      Davide Cenzato and
                      {\relax Zs}uzsanna Lipt{\'{a}}k and
                      Giuseppe Romana},
      title        = {On the Number of Equal-Letter Runs of the Bijective {Burrows-Wheeler}
                      Transform},
      booktitle    = {Proceedings of the 24th Italian Conference on Theoretical Computer
                      Science, Palermo, Italy, September 13-15, 2023},
      volume       = {3587},
      pages        = {129--142},
      year         = {2023}
    }
