% TEST_TSWHIST - Test script for tswHist function validation and performance
%   Validates the correctness of tswHist implementations against exhaustive
%   histcounts computation and measures performance across different variants.
%
%   Tests all variants: builtin, custom-ml, custom-mx, and full MEX implementation.
%   Performs timing benchmarks and validates equivalence between implementations.
%
% Example:
%   run test_tswHist
%
% Other m-files required: tswHist.m, tswHist_mx (MEX), hist_int_mx (MEX)
% Subfunctions: none
% MAT-files required: none
%
% See also: tswHist.m, tswHist_mx.c, tswHist_mx.h, hist_int_mx.c, histcounts
%
% Project: tswHist (https://github.com/cyber-g/tswHist)
%
% License: GNU General Public License v3.0
%
% Author: Germain PHAM
% C2S, Télécom ParisTech, IP Paris
% August 2025; Last revision:

%------------- BEGIN CODE --------------

addpath('..')

% Check if tswHist_mx.mexa64 and hist_int_mx.mexa64 exist
if ~exist('tswHist_mx.mexa64', 'file') || ~exist('hist_int_mx.mexa64', 'file')
    % check if MEX files are compiled at upper level
    if ~exist('../tswHist_mx.mexa64', 'file') || ~exist('../hist_int_mx.mexa64', 'file')
        system('cd .. && make');
    end
end



% Test for tswHist function
n_bins  = 100;  % Number of histogram bins
win_len = 5000; % Length of the sliding window
stride  = 10;  % Stride for the sliding window

x = randn(1, 100000); % Long Gaussian random vector

[histMat_bt, windows_loci_bt, edges_bt] = tswHist(x, n_bins, win_len, stride);
timeit(@() tswHist(x, n_bins, win_len, stride))

[histMat_custml, windows_loci_custml, edges_custml] = tswHist(x, n_bins, win_len, stride,'custom-ml');
timeit(@() tswHist(x, n_bins, win_len, stride,'custom-ml'))

[histMat_custmx, windows_loci_custmx, edges_custmx] = tswHist(x, n_bins, win_len, stride,'custom-mx');
timeit(@() tswHist(x, n_bins, win_len, stride,'custom-mx'))

[histMat_fullmx, windows_loci_fullmx, edges_fullmx] = tswHist_mx(x, n_bins, win_len, stride);
timeit(@() tswHist_mx(x, n_bins, win_len, stride))

histMat_ref = zeros(n_bins, floor((length(x) - win_len + 1) / stride));

% Exhaustive computation for each window
tic
histcounts_edges = (0:n_bins) / n_bins * (max(x) - min(x)) + min(x);
for i = 1:length(windows_loci_bt)
    idx = windows_loci_bt(i):(windows_loci_bt(i)+win_len-1);
    histMat_ref(:, i) = histcounts(x(idx), histcounts_edges);
end
toc

assert(isequal(histMat_bt, histMat_ref), 'Sliding window histograms do not match exhaustive computation.');
assert(isequal(histMat_custml, histMat_ref), 'Custom ML sliding window histograms do not match exhaustive computation.');
assert(isequal(histMat_custmx, histMat_ref), 'Custom MX sliding window histograms do not match exhaustive computation.');
assert(isequal(histMat_fullmx, histMat_ref), 'Full MX sliding window histograms do not match exhaustive computation.');

assert(isequal(windows_loci_bt, windows_loci_custml), 'Window loci do not match between built-in and custom ML.');
assert(isequal(windows_loci_bt, windows_loci_custmx), 'Window loci do not match between built-in and custom MX.');
assert(isequal(windows_loci_bt, windows_loci_fullmx), 'Window loci do not match between built-in and full MX.');

assert(isequal(edges_bt,histcounts_edges), 'Edges do not match between built-in and exhaustive computation.');
assert(isequal(edges_custml, histcounts_edges), 'Edges do not match between custom ML and exhaustive computation.');
assert(isequal(edges_custmx, histcounts_edges), 'Edges do not match between custom MX and exhaustive computation.');
assert(isequal(edges_fullmx, histcounts_edges), 'Edges do not match between full MX and exhaustive computation.');

disp(['All tests in '  mfilename() ' passed successfully!']);
%------------- END OF CODE --------------

