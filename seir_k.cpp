#include <iostream>
#include <sstream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

double b = 1; // Beta
double g = 0.5; // Rate for leaving I state
double v = 0.5; // Rate for leaving E state
const int k = 1000; // Number of E and I states
const long double N = 1000000000;

typedef boost::array< long double , 2+2*k > state_type;

void sir( const state_type &x , state_type &dxdt , double t)
{
	long double I = 0;
	// Compute total I, sum of individual I compartments
	for (int i=k+1; i<2*k+1; i++) {
                I  = I + x[i];
	}
    
	dxdt[0] = -b * x[0] * I/N;
    
	dxdt[1] = b * x[0] * I/N - k*v * x[1]; // First E state
    for (int i=2; i<k+1; i++) { // E states
	    dxdt[i] = k*v * x[i-1] - k*v * x[i];
    }
    
    dxdt[k+1] = k*v * x[k] - k*g * x[k+1]; // First I state
    for (int i=k+2; i < 2*k+1; i++) {
	    dxdt[i] = k*g * x[i-1] - k*g * x[i];
    }
    
    dxdt[2*k+1] = k*g * x[2*k]; // R state
}

void write_sir(const state_type &x , const double t)
{
	string outstr = "";
    ostringstream streamObj;
    streamObj.precision(3);

	streamObj << x[0]/N;
	streamObj << " ";
	
	// The E states:
	long double E=0;
	for (int i=1; i<k+1; i++) { 
		E = E + x[i];
	}
	// The I states:
	long double I = 0;
	for (int i=k+1; i<2*k+1; i++) {
                I  = I + x[i];
        }
	streamObj << E/N;
	streamObj << " ";
	streamObj << I/N;
	streamObj << " ";
	
	streamObj << x[2*k+1]/N;
    outstr = streamObj.str();

	cout << t << " " << outstr << endl;
}

int main(int argc, char **argv)
{
    bool even_start;
    even_start = true; // Evenly divided initially, or all in first I-state?
	if (argc > 2) { // Command line parameters are T and c
		g = 1.0/stof(argv[1]);
		v = (1.0/stof(argv[2])) * g;
		cout << "# Parameters: beta v g k" << endl;
		cout << "# Parameters: " << b << " " << v << " " << g << " " << k << endl;
		state_type x;
		
		x[0]=N-1;
		
		
		for (int i=1; i < x.size(); i++) {
			x[i]=0.0;
		}
		
		if (even_start) {
		    for (int i=1; i < 2*k+1; i++) {
			    x[i]=(N-x[0])/(2*k);
		    }
		}
		else {
		    x[k+1]=N-x[0]; // The first I state!
		}
		
		integrate( sir , x , 0.0 , 3.0*(10+10*10) , 0.01/k , write_sir );

	}
	else {
		cout << "Not enough parameters specified!!" << endl;
		cout << "Please specify T and c as command line parameters, e.g." << endl;
		cout << "sir 5 1.5" << endl;
	}
}

