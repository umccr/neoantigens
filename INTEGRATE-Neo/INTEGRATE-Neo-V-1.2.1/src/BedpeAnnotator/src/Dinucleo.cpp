#include "Dinucleo.h"

int Dinucleo::loadFromFile(char * filename)
{
  spliceVec.clear();
  string line;
  ifstream myfile (filename);
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      splice_t st;
      std::vector<std::string> tmp = my_split(line, '\t');
      st.FiveOne=(tmp[0].c_str())[0];
      st.FiveTwo=(tmp[0].c_str())[1];
      st.ThreeOne=(tmp[1].c_str())[0];
      st.ThreeTwo=(tmp[1].c_str())[1];
     
      spliceVec.push_back(st);      

    }
    myfile.close();
  }
  return 0;
}

int Dinucleo::isCanoSplice(char fiveone, char fivetwo, char threeone, char threetwo)
{
    int res=0;
    for(int i=0;i<this->size();i++)
    {
        splice_t st_here=getSplice(i);
        if(st_here.FiveOne==fiveone && st_here.FiveTwo==fivetwo && st_here.ThreeOne==threeone && st_here.ThreeTwo==threetwo)
         {
             return 1;
         }        
    }
    return res;
}

int Dinucleo::isCanoSplice(string chr5, string strand5, uint32_t junc5, string chr3, string strand3, uint32_t junc3, Reference & ref)
{
    char fiveone='X';
    char fivetwo='X';
    char threeone='X';
    char threetwo='X';

    if(strand5.compare("+")==0)
    {
        fiveone=ref.getRefChar(ref.to_ref_pos(chr5,junc5+1));//+1 and +2 need to check boundary
        fivetwo=ref.getRefChar(ref.to_ref_pos(chr5,junc5+2));
    }
    else
    {
        fiveone=getCharComp(ref.getRefChar(ref.to_ref_pos(chr5,junc5-1)));//+1 and +2 need to check boundary
        fivetwo=getCharComp(ref.getRefChar(ref.to_ref_pos(chr5,junc5-2)));
    }
    if(strand3.compare("+")==0)
    {
        threeone=ref.getRefChar(ref.to_ref_pos(chr3,junc3-2));//+1 and +2 need to check boundary
        threetwo=ref.getRefChar(ref.to_ref_pos(chr3,junc3-1));
    }
    else
    {
        threeone=getCharComp(ref.getRefChar(ref.to_ref_pos(chr3,junc3+2)));//+1 and +2 need to check boundary
        threetwo=getCharComp(ref.getRefChar(ref.to_ref_pos(chr3,junc3+1)));
    }
    return isCanoSplice(fiveone, fivetwo, threeone, threetwo);
}


