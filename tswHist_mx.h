/*
 * tswHist_mx.h - Core routines for fast sliding window histogram computation
 *
 *   Provides helper functions for efficient sliding window histogram calculation
 *   using integer binning and differential updates. Used by the tswHist_mx MEX gateway.
 *
 *   The pushHist and popHist functions incrementally update histogram vectors.
 *   The tswHistSlidingWindow function implements the main sliding window logic.
 *
 *   The core logic is adapted from the essential version of hist_int in:
 *   https://github.com/cyber-g/FastHist
 *
 *   Project: tswHist (https://github.com/cyber-g/tswHist)
 *   License: GNU General Public License v3.0
 *
 *   Author: Germain PHAM
 *   C2S, Télécom ParisTech, IP Paris
 *   August 2025; Last revision:
 */


#include "mex.h"

void pushHist(double *hist_vec, const double *input_int, mwSize len, mwSize n_bins) {
    for (mwSize i = 0; i < len; ++i) {
        int bin = (int)input_int[i];
        if (bin >= 0 && bin < (int)n_bins) {
            hist_vec[bin] += 1;
        }
    }
}

void popHist(double *hist_vec, const double *input_int, mwSize len, mwSize n_bins) {
    for (mwSize i = 0; i < len; ++i) {
        int bin = (int)input_int[i];
        if (bin >= 0 && bin < (int)n_bins) {
            hist_vec[bin] -= 1;
        }
    }
}

void tswHistSlidingWindow(
    double *histMat,
    double *bufferHist,
    const double *input_int,
    const double *strided_windows_loci, 
    mwSize num_windows,
    mwSize win_len,
    mwSize n_bins,
    mwSize stride,
    const double *offsets,
    mwSize input_len
) {
    for (mwSize w = 1; w < num_windows; ++w) {
        // pop indices
        mwSize base_pop = (mwSize)strided_windows_loci[w] - 2; // -1 for 0-based, -1 for previous window
        for (mwSize j = 0; j < stride; ++j) {
            mwSize idx = base_pop + (mwSize)offsets[j];
            if (idx < input_len)
                popHist(bufferHist, &input_int[idx], 1, n_bins);
        }
        // push indices
        mwSize base_push = (mwSize)strided_windows_loci[w] + win_len - 2;
        for (mwSize j = 0; j < stride; ++j) {
            mwSize idx = base_push + (mwSize)offsets[j];
            if (idx < input_len)
                pushHist(bufferHist, &input_int[idx], 1, n_bins);
        }
        // Store
        for (mwSize b = 0; b < n_bins; ++b)
            histMat[b + w * n_bins] = bufferHist[b];
    }
}