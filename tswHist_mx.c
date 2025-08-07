/*
 * tswHist_mx.c - Fast sliding window histogram computation (MEX gateway)
 *
 *   Twin MEX function for tswHist.m
 * 
 *   Implements efficient sliding window histogram calculation for 1D signals using
 *   integer binning and differential updates.
 *
 *   Please read tswHist.m for more information.
 *
 *   Usage from matlab:
 *     [histMat, strided_windows_loci, edges] = tswHist_mx(input, n_bins, win_len, stride)
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
#include "tswHist_mx.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Argument parsing and validation
    if (nrhs < 3 || nrhs > 4)
        mexErrMsgIdAndTxt("tswHist_mx:invalidNumInputs", "Usage: [histMat, strided_windows_loci, edges] = tswHist_mx(input, n_bins, win_len, stride)");

    // Input
    const mxArray *input_mx = prhs[0];
    mwSize n_bins  = (mwSize)mxGetScalar(prhs[1]);
    mwSize win_len = (mwSize)mxGetScalar(prhs[2]);
    // Optional stride, set to 1 if not provided
    mwSize stride  = (nrhs >= 4) ? (mwSize)mxGetScalar(prhs[3]) : 1;

    if (!mxIsDouble(input_mx) || mxIsComplex(input_mx))
        mexErrMsgIdAndTxt("tswHist_mx:inputNotReal", "Input must be a real double vector.");
    mwSize input_len    = mxGetNumberOfElements(input_mx);
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble *input = mxGetDoubles(input_mx);
    #else
        double *input   = mxGetPr(input_mx);
    #endif

    if (n_bins <= 2)
        mexErrMsgIdAndTxt("tswHist_mx:badBins", "Number of bins must be > 2.");
    if (stride >= win_len)
        mexErrMsgIdAndTxt("tswHist_mx:strideWin", "Stride must be less than window length.");

    // Compute strided windows loci
    mwSize num_windows = (input_len - win_len) / stride + 1;
    plhs[1] = mxCreateDoubleMatrix(1, num_windows, mxREAL);
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble *strided_windows_loci = mxGetDoubles(plhs[1]);
    #else
        double *strided_windows_loci   = mxGetPr(plhs[1]);
    #endif
    for (mwSize i = 0; i < num_windows; ++i)
        strided_windows_loci[i] = 1 + i * stride; // MATLAB 1-based

    // Compute histogram bin edges
    double min_val = input[0], max_val = input[0];
    for (mwSize i = 1; i < input_len; ++i) {
        if (input[i] < min_val) min_val = input[i];
        if (input[i] > max_val) max_val = input[i];
    }
    plhs[2] = mxCreateDoubleMatrix(1, n_bins + 1, mxREAL);
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble *edges = mxGetDoubles(plhs[2]);
    #else
        double *edges   = mxGetPr(plhs[2]);
    #endif
    for (mwSize i = 0; i <= n_bins; ++i)
        edges[i] = min_val + (max_val - min_val) * ((double)i / n_bins);

    // Normalize input to integer bins
    double *input_int = (double *)mxMalloc(input_len * sizeof(double));
    for (mwSize i = 0; i < input_len; ++i) {
        double norm = (input[i] - min_val) / (max_val - min_val);
        int bin     = (int)floor(norm * n_bins);
        if (bin == (int)n_bins) bin = n_bins - 1; // Patch for max value
        input_int[i] = (double)bin;
    }

    // Output: histMat
    plhs[0] = mxCreateDoubleMatrix(n_bins, num_windows, mxREAL);
    #if MX_HAS_INTERLEAVED_COMPLEX
        mxDouble *histMat = mxGetDoubles(plhs[0]);
    #else
        double *histMat   = mxGetPr(plhs[0]);
    #endif

    // Compute histogram for the first window
    double *bufferHist = (double *)mxCalloc(n_bins, sizeof(double));
    pushHist(bufferHist, input_int, win_len, n_bins);
    for (mwSize b = 0; b < n_bins; ++b)
        histMat[b] = bufferHist[b];

    // Prepare offsets for pop/push
    double *offsets = (double *)mxMalloc(stride * sizeof(double));
    for (mwSize i = 0; i < stride; ++i)
        offsets[i] = -(double)(stride - 1 - i);

    // Sliding window
    tswHistSlidingWindow(
        histMat,
        bufferHist,
        input_int,
        strided_windows_loci,
        num_windows,
        win_len,
        n_bins,
        stride,
        offsets,
        input_len
    );

    mxFree(input_int);
    mxFree(bufferHist);
    mxFree(offsets);
}