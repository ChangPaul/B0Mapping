/*
 * spha.cpp
 *
 * MEX-file for MATLAB.
 *
 */

#include "mex.h"
#include "spha.h"
#include "mex_err.h"
#include <cmath>

using namespace std;

/*
 * FUNCTION:	spha
 * DESCRIPTION:	Calculate the values of the spherical harmonic
 *				functions from the given coordinate positions.
 * INPUTS:		dims  - dimensions of the output matrix [PxIxN]
 *                      [(no. of pos sets) x (no. of spher harm) x (no. of pos)]
 *              coord - Nx3xP array of coordinate positions.
 *				ind   - vector of spher harmonic indices to be calculated.
 *              normalisation - normalisation flag.
 * OUTPUTS:		output - spher harmonic values at given positions.
 */
void spha (double *output, mwSignedIndex *dims, double *coord, double *ind, bool normalisation)
{
	double val, dist;
	int x;
	
	// Iterate through each set of positions
	for (int p = 0; p != dims[2]; ++p)
	{
		// Iterate through each requested sph harm function
		for (int i = 0; i != dims[1]; ++i)
		{
			int n = dims[0];
			
			// Iterate through all positions
			for (int r = 0; r != n; ++r)
			{
				x = r + (3*n)*p;
				
				// Calculate sph fn value
				val = sphFn [(int)ind[i] - 1](coord [x], coord[x+n], coord[x+2*n]);
		
				// Calculate normalisation if necessary
				if (normalisation)
				{
					dist = sqrt (coord [x    ]*coord [i    ] + 
						         coord [x+  n]*coord [x+  n] +
							     coord [x+2*n]*coord [x+2*n]);
					for (int j = 1; j*j < n; ++j)
						val /= dist;
				}
		
				// Write to output
				x = r + n*i + n*dims[1]*p;
				output[x] = val;
			}
		}
	}
}

//==================================================================
// Main Function
//==================================================================
 
void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int idx, normalisation = false;
	double *coord, *output, *ind;
	size_t n, m;
	mwSignedIndex dims[3];
	const mwSize *pDims;

	// Check the arguments
	int err = 0;
	switch (nrhs)
	{
	case 3: // Normalisation
			if (mxGetNumberOfElements (prhs[2]) != 1 ||
			    !mxIsLogical (prhs[2]))                          err = 0x112;
			else normalisation = (int) *mxGetPr (prhs[2]);
	case 2: // Coordinates
			pDims = mxGetDimensions (prhs[1]);
			if (!mxIsDouble (prhs[1]) || mxIsComplex(prhs[1]))   err = 0x112;
			else if (pDims[1] != 3)                              err = 0x113;
			else coord = mxGetPr (prhs[1]);
			dims[0] = pDims[0];
			dims[2] = (mxGetNumberOfDimensions (prhs[1]) > 2) ? pDims[2] : 1;
    case 1: // Index no.
			m = (int) mxGetM (prhs[0]);
			n = (int) mxGetN (prhs[0]);
			if (!mxIsNumeric (prhs[0]))   					     err = 0x112;
			else if (n != 1 && m != 1)                           err = 0x113;
			else ind = mxGetPr (prhs[0]);
			dims[1] = (m > n) ? m : n;
			break;
	case 0: err = 0x100; break;
	default: err = 0x111;
	}
	
	if (nrhs <  2) err = 0x101;
	if (nlhs != 1) err = 0x200;
	error (err, "spha");
	
	// Create output
	plhs[0] = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	output = mxGetPr (plhs[0]);
	
	// Read probe frequencies and calculate optimum positions
	spha (output, dims, coord, ind, (bool)normalisation);
}