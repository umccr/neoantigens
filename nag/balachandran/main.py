'''
Created on Jul 27, 2017

@author: Marta Luksza, mluksza@ias.edu
'''
from Aligner import Aligner
import sys

def main():        
    '''
    command line parameters:
    neofile - text file with neoantigen data (supplementary data)
    alignmentDirectory - folder with precomputed alignments (SI)
    a - midpoint parameter of the logistic function, alignment score threshold
    k - slope parameter of the logistic function
    '''
    neofile=sys.argv[1]
    alignmentDirectory=sys.argv[2]
    
    a=float(sys.argv[3])
    k=float(sys.argv[4])
    
    #Compute MHC amplitudes for all neoantigens
    f=open(neofile)
    lines=f.readlines()
    Ai={}
    data={}
    samples=set()
    for line in lines[1:]:
        [i,sample,_,_,_,_,mtpeptide,_,_,kdwt,kdmt]=line.strip().split()
        i=int(i)
        data[i]=mtpeptide.upper()
        Ai[i]=float(kdwt)/float(kdmt)
        samples.add(sample)
    f.close()

    #Compute TCR-recognition probabilities for all neoantigens
    aligner=Aligner()    
    for sample in samples:
        xmlpath=alignmentDirectory+"/neoantigens_"+sample+"_iedb.xml"
        aligner.readAllBlastAlignments(xmlpath)    
    aligner.computeR(a, k)    

    #Compute neoantigen quality
    nids=list(Ai.keys())
    nids.sort()
    header=["NeoantigenID","MT.Peptide.Form","NeoantigenQuality",
            "NeoantigenAlignment","IEDB_EpitopeAlignment","AlignmentScore","IEDB_Epitope"]
    header="\t".join(header)
    print(header)
    for i in nids:
        A=Ai[i]
        [R,species,alignment]=aligner.getR(i)
        
        neoAlignment=alignment[0]
        epitopeAlignment=alignment[1]
        score=alignment[2]
        
        l=[i, data[i], A*R, neoAlignment, epitopeAlignment, score, species]
        l="\t".join(map(lambda s: str(s),l))
        print(l)

if __name__ == '__main__':
    if len(sys.argv)!=5:
        print("Run as:")
        print("python src/main.py <Neoantigen_file> <Alignment_directory> <a> <k>")
    else:
        main()
        
