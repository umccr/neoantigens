# Don't source it. Follow carefully.

#VEP_DATA=/Users/vsaveliev/bio/genomes/pcgr/data/grch38/.vep/homo_sapiens/98_GRCh38/ # /g/data3/gx8/extras/umccrise_2019_Aug/genomes/pcgr/data/grch38/.vep
VEP_DATA=/g/data3/gx8/extras/umccrise_genomes/hg38/pcgr/data/grch38/.vep
CONDA_BASE=/g/data/gx8/extras/umccrise_2020_Aug/miniconda/envs

#PYENSEMBL_DIR=/Users/vsaveliev/bio/genomes
PYENSEMBL_CACHE_DIR=/g/data3/gx8/extras/umccrise/genomes/hg38
PIZZLY_ENSEMBL_VERSION=95

# Enviornment
mamba env create -p $CONDA_BASE/nag --file environment.yml
export PATH=$CONDA_BASE/nag/bin:$PATH
pip install -e .

########################
### Install VEP [SKIP that: using the on from PCGR]
#if [ ! -d ${VEP_DATA} ] ; then
#    mkdir ${VEP_DATA}
#    vep_install -a cf -s homo_sapiens -y GRCh38 --DESTDIR ${VEP_DATA} --CACHEDIR ${VEP_DATA}/GRCh38
#    vep_install --AUTO p --NO_HTSLIB --NO_UPDATE --PLUGINS Downstream
#fi
## Or, if using PCGR VEP installation, install the Downstream plugin separately:
git clone https://github.com/Ensembl/VEP_plugins.git
mkdir ${VEP_DATA}/Plugins
cp VEP_plugins/Downstream.pm ${VEP_DATA}/Plugins

########################
### Install pVACtools
mamba env create -p $CONDA_BASE/pvac pvacseq
export PATH=$CONDA_BASE/pvac/bin:$PATH
# override the codebase from our fork:
git clone https://github.com/vladsaveliev/pVACtools
pip install -e pVACtools
pvacseq install_vep_plugin ${VEP_DATA}/Plugins

########################
### Download local IEDB
#IEDB_DIR=/Users/vsaveliev/rjn/extras/iedb
IEDB_DIR=/g/data3/gx8/extras/iedb
if [ ! -d $IEDB_DIR ] ; then
    cd $IEDB_DIR
    wget https://downloads.iedb.org/tools/mhci/2.22.3/IEDB_MHC_I-2.22.3.tar.gz
    wget https://downloads.iedb.org/tools/mhcii/2.22.3/IEDB_MHC_II-2.22.3.tar.gz
    tar -zxvf IEDB_MHC_I-*.tar.gz
    # need python2 to run configure:
    conda create -p $CONDA_BASE/py2 python=2.7
    export PATH=$CONDA_BASE/py2/bin:$PATH
    PY2=$(which python2.7)
    #cd ~/bin
    #ln -s /g/data3/gx8/extras/vlad/miniconda/envs/py2/bin/python2.7
    cd mhc_i
    ./configure

    sed -i '/import pkg_resources/d' method/netmhc-4.0-executable/netmhc_4_0_executable/__init__.py
    sed -i '/import pkg_resources/d' method/netmhcpan-2.8-executable/netmhcpan_2_8_executable/__init__.py
    sed -i 's|/usr/bin/env python|/usr/bin/env python2.7|' method/netmhccons-1.1-executable/netmhccons_1_1_executable/bin/pseudofind
    sed -i 's|/usr/bin/env python|/usr/bin/env python2.7|' method/netmhc-3.4-executable/netmhc_3_4_executable/netMHC

    cd -
    tar -xzvf IEDB_MHC_II-2.17.5.tar.gz
    cd mhc_ii
    python configure.py
fi


########################
### pyensembl (to convert pizzly to bedpe)
pip install "gtfparse>=1.1" "pyensembl>=1.7.2"
if [ ! -d ${PYENSEMBL_CACHE_DIR}/pyensembl ] ; then
    # In 2 steps: first on loging node to make it download the files:
    pyensembl install --release ${PIZZLY_ENSEMBL_VERSION} --species human
    # when it starts `Reading GTF from`, go into a worker node and run again.
fi

