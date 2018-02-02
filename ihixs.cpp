/** \file ihixs.cpp
 *
 * Main/Entry of the program.
 */

#include "src/higgs/ggf_manager/ggf_manager.h"
#include "src/tools/user_interface.h"

void print_logo();

int main(int argc, char** argv)
{
    print_logo();
    UserInterface UI;
    try
    {
        UI.ParseInput(argc,argv);
        GgfManager run(UI);
    }
    catch(const char* s)
    {
        cerr << endl << argv[0] << ": " << s << endl;
    }
    catch(const string s)
    {
        cerr << endl << argv[0] << ": " << s << endl;
    }
    catch(exception &e)
    {
      cerr << endl << argv[0] << ": " << e.what() << endl;
    }
    catch(...)
    {
        cerr<<endl<<"Something went wrong but the exception thrown was not recognized"<<endl;
    }

    return 0;
}



void print_logo()
{
    int vmajor = 2;
    int vminor = 0;
    
    cout<<"\n* * * * * * * * * * * * * * * * * * * * * * * * *";
    cout<<"\n*                                               *";
    cout<<"\n*                                               *";
    cout<<"\n*                                               *";
    cout<<"\n*                ihixs "<<vmajor<<"."<<vminor<<"                      *";
    cout<<"\n*                                               *";
    cout<<"\n*                                               *";
    cout<<"\n*                                               *";
    cout<<"\n* * * * * * * * * * * * * * * * * * * * * * * * *";
    cout<<"\n\n";
}
