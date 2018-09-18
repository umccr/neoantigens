# Don't source it. Follow carefully.

ENSEMBL_VERSION="91"
ENSEMBL_DIR=/g/data3/gx8/extras/ensembl
export VEP_DATA=${ENSEMBL_DIR}/vep_data

# Enviornment
conda env create -p nag --file environment.yml
conda activate nag
pip install -e .


########################
### Install VEP
if [ ! -d ${VEP_DATA} ] ; then
    mkdir ${VEP_DATA}
    vep_install -a cf -s homo_sapiens -y GRCh38 --DESTDIR ${VEP_DATA} --CACHEDIR ${VEP_DATA}/GRCh38
    vep_install --AUTO p --NO_HTSLIB --NO_UPDATE --PLUGINS Downstream
fi


########################
### Install pVACtools
pip install pvactools
# override the codebase from our fork:
git clone https://github.com/vladsaveliev/pVACtools ; cd pVACtools ; pip install . ;  cd ..
pvacseq install_vep_plugin $VEP_DATA/Plugins


########################
### Download local IEDB
export IEDB_DIR=/g/data3/gx8/extras/iedb
if [ ! -d $IEDB_DIR ] ; then
    cd $IEDB_DIR
    # Download from [http://tools.iedb.org/mhcii/download](http://tools.iedb.org/mhcii/download/):
    # - click `MHC Class I`
    # - click `To download the tools in tar.gz format: Agree and Download`
    # - do the the same for `MHC Class II`.
    # - Uncompressing and install:
    tar -zxvf IEDB_MHC_I-2.19.1.tar.gz
    # need python2 to run configure:
    conda create -n py2 python=2.7
    conda activate py2
    PY2=$(which python2.7)
    #cd ~/bin
    #ln -s /g/data3/gx8/extras/vlad/miniconda/envs/py2/bin/python2.7
    cd -
    cd mhc_i
    ./configure

    sed -i '/import pkg_resources/d' method/netmhc-4.0-executable/netmhc_4_0_executable/__init__.py
    sed -i '/import pkg_resources/d' method/netmhcpan-2.8-executable/netmhcpan_2_8_executable/__init__.py
    sed -i 's|/usr/bin/env python|/usr/bin/env python2.7|' method/netmhccons-1.1-executable/netmhccons_1_1_executable/bin/pseudofind
    sed -i 's|/usr/bin/env python|/usr/bin/env python2.7|' method/netmhc-3.4-executable/netmhc_3_4_executable/netMHC

    tar -xzvf IEDB_MHC_II-2.17.5.tar.gz
    cd mhc_ii
    ./configure.py

    conda activate pvac
fi


########################
### pyensembl (to convert pizzly to bedpe)
pip install "gtfparse>=1.1" pyensembl
# override the codebase from our fork:
git clone https://github.com/vladsaveliev/pyensembl ; cd pyensembl ; pip install . ;  cd ..
export PYENSEMBL_CACHE_DIR=$ENSEMBL_DIR
if [ ! -d $PYENSEMBL_CACHE_DIR/pyensembl ] ; then
    # In 2 steps: first on loging node to make it download the files:
    pyensembl install --release $ENSEMBL_VERSION --species human
    # when it starts `Reading GTF from`, go into a worker node and run again.
fi


########################
### bam-readcount (for pvacseq coverage of rnaseq)
# bam-readcount cannot install into the same conda env, thus we are creating a separate one:
conda create -n bam-readcount -c bioconda -c conda-forge bam-readcount
conda activate bam-readcount
BAM_READCOUNT=$(which bam-readcount)
conda activate nag
cp $BAM_READCOUNT $CONDA_PREFIX/bin