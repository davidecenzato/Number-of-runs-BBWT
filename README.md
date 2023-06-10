# BBWT_paper

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
