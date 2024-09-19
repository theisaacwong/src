# snakemake --forceall --rulegraph | dot -Tpdf > dag.pdf
#  chmod a-w -R

# Module organization
 /net/eichler/vol26/7200/software/modules-sw
 /net/eichler/vol26/7200/software/modules-repo/CentOS7

# Add to .bash_profile
```bash
export MODSW=/net/eichler/vol26/7200/software/modules-sw
export MODREP=/net/eichler/vol26/7200/software/modules-repo/${MODULES_REL}
umask 002
```

# Methods for installing modules
## Conda Envsssh lige
## https://github.com/deeptools/HiCExplorer
```bash
MODNAME="hicexplorer"
MODVERSION="3.7.3"

mkdir -p ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH}
conda create -p ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH} -c bioconda -c conda-forge hicexplorer
```


## Virtual environment
## https://github.com/open2c/cooltools
```bash
MODNAME="cooltools"
MODVERSION="0.6.1"

mkdir -p ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH}
cd ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH}
module load miniconda/4.12.0
# First try
python3 -m virtualenv --system-site-packages ct_env
source ct_env/bin/activate
pip install cooltools==${MODVERSION}

# Second try
python3 -m virtualenv ct_env
source ct_env/bin/activate
pip install cooltools==${MODVERSION}
```


## pip install from git
## https://github.com/vaquerizaslab/fanc
```bash
MODNAME="fan-c"
MODVERSION="202312-55ec737"

mkdir -p ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH}
module load hdf5/1.8.13-cxxenabled miniconda/4.12.0
cd ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH}
git clone git@github.com:vaquerizaslab/fanc.git
cd fanc
pip install .
```

## Install from source
```bash
MODNAME="hifiasm"
MODVERSION="0.19.8"

mkdir -p ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH}
cd ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH}
cd hifiasm && git checkout 0.19.8 && make
```


## Pre-built binaries
## https://github.com/nanoporetech/dorado
```bash
MODNAME="dorado"
MODVERSION="0.5.0"

mkdir -p ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH}
cd ${MODSW}/${MODNAME}/${MODVERSION}/${ARCHPATH}
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.5.0-linux-x64.tar.gz
tar zxvf dorado-0.5.0-linux-x64.tar.gz
ln -s dorado-0.5.0-linux-x64/bin/
ln -s dorado-0.5.0-linux-x64/lib/
```


# Module files
```
#%Module1.0

prereq   modules-eichler
conflict miniconda

set      MOD_ROOT     $env(MOD_EICHLERSW)
set      MOD_NAME     pgrtk
set      MOD_VERSION  0.5.1
set      MOD_DIR      $MOD_ROOT/$MOD_NAME/$MOD_VERSION/$env(MODULES_OS)/$env(MODULES_REL)/$env(MODULES_MACH)

module-whatis "$MOD_NAME $MOD_VERSION."

if ![file isdirectory $MOD_DIR] {
    puts stderr "Couldn't find directory $MOD_DIR"
    
    return 1
}

prepend-path    PATH    $MOD_DIR/bin/
#prepend-path    LD_LIBRARY_PATH    $MOD_DIR/dist/lib

setenv   MOD_MINICONDA_BASE $MOD_DIR
```
