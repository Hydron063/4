#include <iostream>
#include "minimisation.h"

using namespace std;

int main()
{
    int pas=1,n=2;  //n - nombre de variables de la fonction
    vector <double> x(n),s(n);
    for (int i=0;i<n;++i)
        cin>>x[i];
    demontrer(minMultidim(x,0.001,0.001));
}
