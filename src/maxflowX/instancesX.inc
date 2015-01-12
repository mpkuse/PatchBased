#include "graphX.h"

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

// Instantiations: <captype, tcaptype, flowtype>
// IMPORTANT: 
//    flowtype should be 'larger' than tcaptype 
//    tcaptype should be 'larger' than captype

template class GraphX<int,int,int>;
//template class GraphX<short,int,int>;
//template class GraphX<float,float,float>;
//template class GraphX<double,double,double>;

