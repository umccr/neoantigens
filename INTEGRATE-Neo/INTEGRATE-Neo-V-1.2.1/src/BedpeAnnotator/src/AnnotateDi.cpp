#include "AnnotateDi.h"

int AnnotateDi::getCanonicalDi(bedpe_t & bp, Dinucleo & di, Reference & ref)
{
    string chr5;
    string strand5;
    uint32_t junc5;
    string chr3;
    string strand3;
    uint32_t junc3;

    chr5=bp.chr1;
    strand5=bp.strand1;
    if(strand5.compare("+")==0)
    {
        junc5=bp.end1;
    }
    else
    {
        junc5=bp.start1+1;
    }

    chr3=bp.chr2;
    strand3=bp.strand2;
    if(strand3.compare("+")==0)
    {
        junc3=bp.start2+1;
    }
    else
    {
        junc3=bp.end2;
    }

    return di.isCanoSplice(chr5, strand5, junc5, chr3, strand3, junc3, ref);

}

int AnnotateDi::annotateDi(Bedpe & bedpe, Dinucleo & di, Reference & ref)
{

    for(int i=0;i<bedpe.size();i++)
    {
        bedpe_t bpe=bedpe.getBedpe(i);
        bedpe.insertCol(i,to_string(static_cast<long long>(getCanonicalDi(bpe,di,ref))));
    }

    return 0;

}
