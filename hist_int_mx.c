/*
 * hist_int_mx.c - Efficient integer histogram for MATLAB (MEX gateway)
 *
 *   Implements a fast histogram for integer-valued input, used in sliding window
 *   histogram computations. This MEX function is an exact C transcription of the
 *   local hist_int MATLAB function, which is itself adapted from the core logic of:
 *   https://github.com/cyber-g/FastHist (only the essential logic is retained).
 *
 *   Usage:
 *     hist_vec = hist_int_mx(input_int, n_bins)
 *
 *   Project: tswHist (https://github.com/cyber-g/tswHist)
 *   License: GNU General Public License v3.0
 *
 *   Author: Germain PHAM
 *   C2S, Télécom ParisTech, IP Paris
 *   August 2025; Last revision:
 */

#include "mex.h"
#include "matrix.h"
#include "tswHist_mx.h"

/* Gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /*
     * Usage:
     *   hist_vec = hist_int_mx(input_int, n_bins)
     */

    if (nrhs == 2) {
        // hist_int_mx(input_int, n_bins)
        // input_int: vector, n_bins: scalar
        if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
            mexErrMsgIdAndTxt("hist_int_mx:inputNotReal", "input_int must be a real double vector.");
        if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) || mxGetNumberOfElements(prhs[1]) != 1)
            mexErrMsgIdAndTxt("hist_int_mx:binsNotScalar", "n_bins must be a real scalar.");

        mwSize n_bins = (mwSize)mxGetScalar(prhs[1]);
        mwSize len = mxGetNumberOfElements(prhs[0]);
        #if MX_HAS_INTERLEAVED_COMPLEX
            mxDouble *input_int = mxGetDoubles(prhs[0]);
        #else
            double *input_int = mxGetPr(prhs[0]);
        #endif

        // hist_int
        plhs[0] = mxCreateDoubleMatrix(1, n_bins, mxREAL); // array is initialized to zero by mxCreateDoubleMatrix
        #if MX_HAS_INTERLEAVED_COMPLEX
            mxDouble *hist_vec = mxGetDoubles(plhs[0]);
        #else
            double *hist_vec = mxGetPr(plhs[0]);
        #endif
        pushHist(hist_vec, input_int, len, n_bins);

    } else {
        mexErrMsgIdAndTxt("hist_int_mx:invalidNumInputs",
            "Usage: hist_vec = hist_int_mx(input_int, n_bins)");
    }
}