#This example is using GRCh38.
#The GRCh38.fa is not provided in this directory due to the size of it.  
#If you install INTEGRATE-Neo and have the setup.ini set correctly, you should be able to run this:
#(but make sure to change the "/PATH/TO"s to proper paths before running).

python /PATH/TO/integrate-neo.py -1 rd1.fq -2 rd2.fq -f fusions.bedpe -r /PATH/TO/GRCh38.fa -g /PATH/TO/GRCh38.genePhred -k

#and get the output result.bedpe in the fusion_antigen_out directory.
