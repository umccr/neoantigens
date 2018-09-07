#include "Bedpe.h"


int Bedpe::loadFromFile(char * filename)
{
  bedpevec.clear();
  string line;
  ifstream myfile (filename);
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
        bedpe_t bt;
        std::vector<std::string> tmp = my_split(line, '\t');
        bt.chr1=tmp[0];
        bt.start1=atol(tmp[1].c_str());
        bt.end1=atol(tmp[2].c_str());
        bt.chr2=tmp[3];
        bt.start2=atol(tmp[4].c_str());
        bt.end2=atol(tmp[5].c_str());
        bt.name=tmp[6];
        bt.score=atof(tmp[7].c_str()); 
        bt.strand1=tmp[8][0];
        bt.strand2=tmp[9][0];
        if(tmp.size()>10)
        {
            vector<string>::iterator it;
            for(it=tmp.begin()+10;it!=tmp.end();it++)
                bt.others.push_back(*it);
        }
        bedpevec.push_back(bt);
    }
    myfile.close();
  }    

  return 0;
}

int Bedpe::getPos(bedpe_t & tt, uint32_t & tt_pos5, uint32_t & tt_pos3)
{
    if(tt.strand1.compare("+")==0)
    {
        tt_pos5=tt.end1;
    }
    else
    {
        tt_pos5=tt.start1+1;
    }

    if(tt.strand2.compare("+")==0)
    {
        tt_pos3=tt.start2+1;
    }
    else
    {
        tt_pos3=tt.end2;
    }

    return 0;
}


int Bedpe::printBedpe(char * file)
{
    
    std::ofstream ofs;
    ofs.open (file, std::ofstream::out);
    
    for(int i=0;i<bedpevec.size();i++)
    {
        bedpe_t bt=bedpevec[i];    
        ofs<<bt.chr1<<"\t";
        if(bt.start1==-1)
            ofs<<"-1"<<"\t";
        else
            ofs<<bt.start1<<"\t";
        if(bt.end1==-1)
            ofs<<"-1"<<"\t";
        else
            ofs<<bt.end1<<"\t";
        ofs<<bt.chr2<<"\t";
        if(bt.start2==-1)
            ofs<<"-1"<<"\t";
        else          
            ofs<<bt.start2<<"\t";
        if(bt.end2==-1)
            ofs<<"-1"<<"\t";
        else
            ofs<<bt.end2<<"\t";
        ofs<<bt.name<<"\t";
        ofs<<bt.score<<"\t";
        ofs<<bt.strand1<<"\t";
        ofs<<bt.strand2;
        if(bt.others.size()==0)
           ofs<<"\n";
        else 
        {
           ofs<<"\t";
           for(int j=0;j<bt.others.size()-1;j++)        
           {
               if((bt.others[j]).compare("")==0)
                   ofs<<"NA"<<"\t";
               else
                   ofs<<bt.others[j]<<"\t";
           }
               if((bt.others[bt.others.size()-1]).compare("")==0)
                   ofs<<"NA"<<endl;
               else
                   ofs<<bt.others[bt.others.size()-1]<<endl;
        }
    } 
    ofs.close();

    return 0;
}

bool my_sort_bedpe(const bedpe_t & bt1, const bedpe_t & bt2)
{
    string chr1_1=bt1.chr1;
    string chr1_2=bt1.chr2;
    string str1_1=bt1.strand1;
    string str1_2=bt1.strand2;

    string chr2_1=bt2.chr1;
    string chr2_2=bt2.chr2;
    string str2_1=bt2.strand1;
    string str2_2=bt2.strand2;    

    uint32_t bt1_pos5,bt1_pos3;
    uint32_t bt2_pos5,bt2_pos3;
    
    if(bt1.strand1.compare("+")==0)
    {
        bt1_pos5=bt1.end1;
    }
    else
    {
        bt1_pos5=bt1.start1+1;
    }

    if(bt1.strand2.compare("+")==0)
    {
        bt1_pos3=bt1.start2+1;
    }
    else
    {
        bt1_pos3=bt1.end2;
    }

    if(bt2.strand1.compare("+")==0)
    {
        bt2_pos5=bt2.end1;
    }
    else
    {
        bt2_pos5=bt2.start1+1;
    }

    if(bt2.strand2.compare("+")==0)
    {
        bt2_pos3=bt2.start2+1;
    }
    else
    {
        bt2_pos3=bt2.end2;
    }


    if(chr1_1.compare(chr2_1)<0)
    {
        return true;
    }
    else if(chr1_1.compare(chr2_1)==0)
    {
        if(chr1_2.compare(chr2_2)<0)
        {
            return true;
        }
        else if(chr1_2.compare(chr2_2)==0)
        {
            if(str1_1.compare(str2_1)<0)
            {
                return true;
            }
            else if(str1_1.compare(str2_1)==0)
            {
                if(str1_2.compare(str2_2)<0)
                {
                    return true;
                }
                else if(str1_2.compare(str2_2)==0)
                {
                    if(bt1_pos5<bt2_pos5)
                    {
                        return true;
                    }
                    else if(bt1_pos5==bt2_pos5)
                    {
                        if(bt1_pos3<bt2_pos3)
                            return true;
                        else 
                            return false;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    return false;
                }
            }
            else 
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

int Bedpe::uniq()
{
    sort(bedpevec.begin(),bedpevec.end(),my_sort_bedpe);
    vector<bedpe_t> bedpevec2;
    for(int i=0;i<=bedpevec.size();i++)
    {
        if(i==0)
        {
            bedpevec2.push_back(bedpevec[0]);
        }
        if(i>=1)
        {
            bool add=true;
            if (bedpevec[i].chr1.compare(bedpevec[i-1].chr1)==0 && bedpevec[i].chr2.compare(bedpevec[i-1].chr2)==0 && bedpevec[i].strand1.compare(bedpevec[i-1].strand1)==0 && bedpevec[i].strand2.compare(bedpevec[i-1].strand2)==0)
            {
                int delta1=0;
                int delta2=0;

                uint32_t bt1_pos5,bt1_pos3;
                uint32_t bt2_pos5,bt2_pos3;

                if(bedpevec[i-1].strand1.compare("+")==0)
                {
                    bt1_pos5=bedpevec[i-1].end1;
                }
                else
                {
                    bt1_pos5=bedpevec[i-1].start1+1;
                }

                if(bedpevec[i-1].strand2.compare("+")==0)
                {
                    bt1_pos3=bedpevec[i-1].start2+1;
                }
                else
                {
                    bt1_pos3=bedpevec[i-1].end2;
                }

                if(bedpevec[i].strand1.compare("+")==0)
                {
                    bt2_pos5=bedpevec[i].end1;
                }
                else
                {
                    bt2_pos5=bedpevec[i].start1+1;
                }

                if(bedpevec[i].strand2.compare("+")==0)
                {
                    bt2_pos3=bedpevec[i].start2+1;
                }
                else
                {
                    bt2_pos3=bedpevec[i].end2;
                }

                if (bedpevec[i].strand1.compare("+")==0)
                    delta1 = bt2_pos5-bt1_pos5;
                else
                    delta1 = bt1_pos5-bt2_pos5;
                if (bedpevec[i].strand2.compare("+")==0)
                    delta2 = bt2_pos3-bt1_pos3;
                else
                    delta2 = bt1_pos3-bt2_pos3;

                if(delta1 == delta2)
                {
                    add=false;
                }
            }
            
            if(add==true)
            {
                bedpevec2.push_back(bedpevec[i]);
            }
        }
    }
    bedpevec.clear();
    bedpevec=bedpevec2;
    return 0;
}

int Bedpe::insertCol(int index, string str)
{
    if(index>=bedpevec.size())
    {
        cerr<<index<<" > bedpe rows"<<endl;
        exit(1);
    }
    bedpevec[index].others.push_back(str);
    return 0;    
}

int Bedpe::isBedpeSame(bedpe_t & a, bedpe_t & b)
{
    string chr5pA,chr3pA;
    uint32_t bt1_pos5,bt1_pos3;
    string strand5pA,strand3pA;
    string chr5pB,chr3pB;
    uint32_t bt2_pos5,bt2_pos3;
    string strand5pB,strand3pB;
    

    chr5pA=a.chr1;
    chr3pA=a.chr2;
    
    chr5pB=b.chr1;
    chr3pB=b.chr2;

    getPos(a,bt1_pos5,bt1_pos3);
    getPos(b,bt2_pos5,bt2_pos3);


    strand5pA=a.strand1;
    strand3pA=a.strand2;
    
    strand5pB=b.strand1;
    strand3pB=b.strand2;

    int same=0;
    if (chr5pA.compare(chr5pB)==0 && chr3pA.compare(chr3pB)==0 && strand5pA.compare(strand5pB)==0 && strand3pA.compare(strand3pB)==0)
    {
        int delta1=0;
        int delta2=0;

        if (a.strand1.compare("+")==0)
            delta1 = bt2_pos5-bt1_pos5;
        else
            delta1 = bt1_pos5-bt2_pos5;
        if (a.strand2.compare("+")==0)
            delta2 = bt2_pos3-bt1_pos3;
        else
            delta2 = bt1_pos3-bt2_pos3;

        if(delta1 == delta2)
        {
            same=1;
        }
    }
    return same;
}
