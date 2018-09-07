#include "CutterByRule.h"

int CutterByRule::cut(Bedpe & in, Bedpe & out, string rule, string file)
{
    string from = file+".from";
    char * fromfile=(char *)from.c_str();
    in.printBedpe(fromfile);
    string awkStr="awk '"+rule+"{print}' "+from+" > "+file;
    system(awkStr.c_str());
    
    out.loadFromFile((char *)(file.c_str()));
    system(("rm "+from).c_str());  

    return 0;
}
