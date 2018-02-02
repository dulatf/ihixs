#include "iostream"
#include "iomanip"
using namespace std;
#include "result_pair.h"


ostream& operator<<(ostream& stream, const ResultPair& rs)
{
    if (fabs(rs.val())>1e-3)
    {
        stream<<setw(11)<<setprecision(5)<<right<<fixed<<rs.val()
            <<scientific<<setprecision(0)<<left<<"["<<rs.err()<<"] "
            <<setprecision(16)<<fixed;
    }
    else
    {
        stream<<setw(11)<<setprecision(5)<<right<<scientific<<rs.val()
        <<scientific<<setprecision(0)<<left<<"["<<rs.err()<<"] "
        <<setprecision(16)<<fixed;
    }
    return stream;
}