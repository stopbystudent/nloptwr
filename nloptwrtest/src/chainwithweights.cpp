
#include "nloptwrtest/chainwithweights.h"

#include "nloptwr/nloptwr.h"
#include "nloptwr/optfktnbase.h"
#include "nloptwr/optfktnclass.h"

#include <cmath>
#include <vector>

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

namespace opttest
{
    
    // static 
    const double ChainWithWeights::PI = std::atan(1.0)*4.0;
    
// static 
    const double ChainWithWeights::epsilon=1E-8;

ChainWithWeights::ChainWithWeights( int dim, double xN, double hN, double lM, double lI, const std::vector<double>& myWeights)
  :
  oif::OptFknClass(),
  nDim ( dim ),
  x0(0.0),
  h0(0.0),
  xn(xN),
  hn(hN),
  lm(lM),
  li(lI),
  weights(dim+2, 0.0),
  xc(dim+2, 0.0),
  hc(dim+2, 0.0)
{
    for (size_t i=0; (i<dim)&&(i<myWeights.size()); i++) {
        weights[i+1] = myWeights.at(i);
    }
    initialize ();
    validate();
}

// virtual
void ChainWithWeights::initialize () //  double lb, double ub, double xInit )
{
    init( nDim, mDim, -0.5*PI, 0.5*PI, 0.0 );

    xc[0]=x0;
    hc[0]=h0;
    for (size_t i=1; i<=(nDim+1); i++) {
        xc[i] = xc[i-1];
        hc[i] = hc[i-1];
    }    
    xc[nDim+1]=xn;
    hc[nDim+1]=hn;
    
    for (size_t i=0; i<nDim; i++) {
        setXInit(i, 0.0);
    }
}

// virtual
ChainWithWeights::~ChainWithWeights() {}

// The clone function is essential for parallelization
oif::OptFknBase* ChainWithWeights::clone() const
{
  return ( new ChainWithWeights ( *this ) );
}

// virtual
double ChainWithWeights::optFktn ( const std::vector<double>& x, std::vector<double>& c )
{
    // lefe side
    hc[0] = h0;
    xc[0] = x0;

    xc[nDim+1] = xn;
    hc[nDim+1] = hn;

    c[0] = 0.0; 
    c[1] = 0.0; 
    
    // copy variable x[i] into hc[i]
    for (size_t i=0; i<nDim; i++) {
        xc[i+1] = xc[i] + li*cos(x[i]);
        hc[i+1] = hc[i] + li*sin(x[i]);
    }

    double diffL = ((xc[nDim] - xn)*(xc[nDim] - xn) + (hc[nDim] - hn)*(hc[nDim] - hn)) - li*li;
    
    // equality constraint as two inequality constraints
    c[0] = (diffL > epsilon )?  1000.0*diffL : 0.0; 
    c[1] = (diffL < epsilon )? -1000.0*diffL : 0.0; 

    // W_pot/g => m*h => min
    double w_g=0.0;
    double w_g_chain=0.0;
    double w_g_weights=0.0;
    
    for (size_t i=1; i<=(nDim+1); i++) {
        w_g_chain += (hc[i-1]+hc[i])*0.5;
        w_g_weights += hc[i]*weights[i];
    }
    w_g_chain *= (lm*li);
    
    w_g = (w_g_chain + w_g_weights);

    /*
    cout << " " << endl;
    for (size_t i=0; i<=(nDim+1); i++) {
        cout << "pos("<< setw(3) << i << ") = { " << setw(7) << xc[i] << ", " << setw(7) << hc[i] << " } " << endl;
    }
    
    cout << "w_g=" << setw(7) << w_g << ", diffL=" << setw(7) << diffL << endl;
    */
    
  return w_g;
}

} // namespace opttest


