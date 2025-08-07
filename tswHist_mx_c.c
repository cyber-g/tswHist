/*
 * tswHist_mx.c - Fast sliding window histogram computation (MEX gateway)
 *
 *   Twin MEX function for tswHist.m
 * 
 *   Implements efficient sliding window histogram calculation for 1D signals using
 *   integer binning and differential updates.
 * 
 *   Compared tswHist_mx.c, it uses a pure C implementation of the sliding
 *   window histogram computation. This done as a basis for projects which do
 *   not use MATLAB's MEX interface. 
 *
 *   Please read tswHist.m for more information.
 *
 *   Usage from matlab:
 *     [histMat, strided_windows_loci, edges] = tswHist_mx_c(input, n_bins, win_len, stride)
 *
 *   Inputs:
 *     input    - Input vector (real double, 1D)
 *     n_bins   - Number of histogram bins (integer > 2)
 *     win_len  - Sliding window length
 *     stride   - Stride for sliding window (default: 1)
 *
 *   Outputs:
 *     histMat              - n_bins x num_windows matrix of histograms
 *     strided_windows_loci - Start indices of each window (1-based)
 *     edges                - Bin edges used for histogramming
 *
 *   See also: tswHist.m, hist_int_mx.c
 *
 *   Project: tswHist (https://github.com/cyber-g/tswHist)
 *
 *   License: GNU General Public License v3.0
 * 
 *   Author: Germain PHAM
 *   C2S, Télécom Paris, IP Paris
 *   August 2025; Last revision:
 */

#include "mex.h"
#include <math.h>
#include "tswHist.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Argument parsing and validation
    if (nrhs < 3 || nrhs > 4)
        mexErrMsgIdAndTxt("tswHist_mx:invalidNumInputs", "Usage: [histMat, strided_windows_loci, edges] = tswHist_mx_c(input, n_bins, win_len, stride)");

    // Input
    const mxArray *input_mx = prhs[0];
    mwSize n_bins  = (mwSize)mxGetScalar(prhs[1]);
    mwSize win_len = (mwSize)mxGetScalar(prhs[2]);
    mwSize stride  = (nrhs >= 4) ? (mwSize)mxGetScalar(prhs[3]) : 1;

    if (!mxIsDouble(input_mx) || mxIsComplex(input_mx))
        mexErrMsgIdAndTxt("tswHist_mx:inputNotReal", "Input must be a real double vector.");
    mwSize input_len = mxGetNumberOfElements(input_mx);
#if MX_HAS_INTERLEAVED_COMPLEX
    const double *input = mxGetDoubles(input_mx);
#else
    const double *input = mxGetPr(input_mx);
#endif

    if (n_bins <= 2)
        mexErrMsgIdAndTxt("tswHist_mx:badBins", "Number of bins must be > 2.");
    if (stride >= win_len)
        mexErrMsgIdAndTxt("tswHist_mx:strideWin", "Stride must be less than window length.");

    // Compute number of windows
    mwSize num_windows = (input_len - win_len) / stride + 1;

    // Allocate outputs
    plhs[0] = mxCreateDoubleMatrix(n_bins, num_windows, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1, num_windows, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, n_bins + 1, mxREAL);

#if MX_HAS_INTERLEAVED_COMPLEX
    double *histMat = mxGetDoubles(plhs[0]);
    double *strided_windows_loci = mxGetDoubles(plhs[1]);
    double *edges = mxGetDoubles(plhs[2]);
#else
    double *histMat = mxGetPr(plhs[0]);
    double *strided_windows_loci = mxGetPr(plhs[1]);
    double *edges = mxGetPr(plhs[2]);
#endif

    // Call pure C implementation
    tswHist(
        input, input_len,
        n_bins, win_len, stride,
        histMat,
        strided_windows_loci,
        edges
    );

}