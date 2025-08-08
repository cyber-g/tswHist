/*
 * tswHist.h - Core routines for fast sliding window histogram computation (pure C version)
 *
 *   Provides helper functions for efficient sliding window histogram calculation
 *   using integer binning and differential updates. This header is pure C and
 *   does not depend on MATLAB or MEX headers.
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

#ifndef TSWHIST_H
#define TSWHIST_H

#include <stddef.h> // for size_t
#include <stdlib.h> // for malloc, free
#include <math.h> // for floor

void pushHist(double *hist_vec, const double *input_int, size_t len, size_t n_bins) {
    for (size_t i = 0; i < len; ++i) {
        int bin = (int)input_int[i];
        if (bin >= 0 && bin < n_bins) {
            hist_vec[bin] += 1;
        }
    }
}

void popHist(double *hist_vec, const double *input_int, size_t len, size_t n_bins) {
    for (size_t i = 0; i < len; ++i) {
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
    size_t num_windows,
    size_t win_len,
    size_t n_bins,
    size_t stride,
    const double *offsets,
    size_t input_len
) {
    for (size_t w = 1; w < num_windows; ++w) {
        // pop indices
        size_t base_pop = (size_t)strided_windows_loci[w] - 2; // -1 for 0-based, -1 for previous window
        for (size_t j = 0; j < stride; ++j) {
            size_t idx = base_pop + (size_t)offsets[j];
            if (idx < input_len)
                popHist(bufferHist, &input_int[idx], 1, n_bins);
        }
        // push indices
        size_t base_push = (size_t)strided_windows_loci[w] + win_len - 2;
        for (size_t j = 0; j < stride; ++j) {
            size_t idx = base_push + (size_t)offsets[j];
            if (idx < input_len)
                pushHist(bufferHist, &input_int[idx], 1, n_bins);
        }
        // Store
        for (size_t b = 0; b < n_bins; ++b)
            histMat[b + w * n_bins] = bufferHist[b];
    }
}

void tswHist(
    const double *input_norm, size_t input_len,
    size_t n_bins, size_t win_len, size_t stride,
    double *histMat,         // [n_bins x num_windows] output
    double *strided_windows_loci, // [num_windows] output
    double *edges            // [n_bins+1] output
) {
    // Compute number of windows
    size_t num_windows = (input_len - win_len) / stride + 1;

    // Compute strided windows loci (maintain 1-based for MATLAB compatibility)
    for (size_t i = 0; i < num_windows; ++i)
        strided_windows_loci[i] = (double)(i * stride + 1); // 1-based

    // Compute the edges for the histogram bins. Bin edges are set to be between
    // 0 and 1 exactly here, 0 and 1 are included in the edges
    for (size_t i = 0; i <= n_bins; ++i)
        edges[i] = 1.0 * ((double)i / n_bins);

    // Normalize input to integer bins
    double *input_int = (double *)calloc(input_len, sizeof(double));
    for (size_t i = 0; i < input_len; ++i) {
        //  The normalization is left outside this function for more flexibility
        //  The input vector is expected to be included in [0,1] (not
        //  necessarily exactly occupying this range)
        int bin = (int)floor(input_norm[i] * n_bins);
        if (bin == (int)n_bins) bin = n_bins - 1; // Patch for max value
        input_int[i] = (double)bin;
    }

    // Compute histogram for the first window
    double *bufferHist = (double *)calloc(n_bins, sizeof(double));
    pushHist(bufferHist, input_int, win_len, n_bins);
    for (size_t b = 0; b < n_bins; ++b)
        histMat[b] = bufferHist[b];

    // Prepare offsets for pop/push
    double *offsets = (double *)calloc(stride, sizeof(double));
    for (size_t i = 0; i < stride; ++i)
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

    free(input_int);
    free(bufferHist);
    free(offsets);
}

#endif // TSWHIST_H