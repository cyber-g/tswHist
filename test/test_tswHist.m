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
% Other m-files required: tswHist.m, tswHist_mx (MEX), hist_int_mx (MEX), tswHist_mx_c (MEX)
% Subfunctions: none
% MAT-files required: none
%
% See also: tswHist.m, tswHist_mx.c, tswHist_mx.h, hist_int_mx.c, tswHist_mx_c.c, histcounts
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

test = 'gaussian'; % Choose 'gaussian' or 'uniform' for the type of input vector to test

% Test for tswHist function
n_bins  = 100;  % Number of histogram bins
win_len = 5000; % Length of the sliding window
stride  = 10;  % Stride for the sliding window

switch test
    case 'gaussian'
        % Generate a long Gaussian random vector normalized to [0, 1]
        x = randn(1, 100000); % Long Gaussian random vector
        x = (x - min(x)) / (max(x) - min(x)); % Normalize to [0, 1]
        % then the data exactly contains at least one 0 and one 1
    case 'uniform'
        % Generate a long uniform random vector normalized to [0, 1]
        x = rand(1, 100000); % Long uniform random vector
        % 0 and 1 are not necessarily in the samples
    otherwise
        error('Unknown test type. Use "gaussian" or "uniform".');
end


[histMat_bt, windows_loci_bt, edges_bt] = tswHist(x, n_bins, win_len, stride);
timeit(@() tswHist(x, n_bins, win_len, stride))

[histMat_custml, windows_loci_custml, edges_custml] = tswHist(x, n_bins, win_len, stride,'custom-ml');
timeit(@() tswHist(x, n_bins, win_len, stride,'custom-ml'))

[histMat_custmx, windows_loci_custmx, edges_custmx] = tswHist(x, n_bins, win_len, stride,'custom-mx');
timeit(@() tswHist(x, n_bins, win_len, stride,'custom-mx'))

[histMat_fullmx, windows_loci_fullmx, edges_fullmx] = tswHist_mx(x, n_bins, win_len, stride);
timeit(@() tswHist_mx(x, n_bins, win_len, stride))

[histMat_mx_c, windows_loci_mx_c, edges_mx_c] = tswHist_mx_c(x, n_bins, win_len, stride);
timeit(@() tswHist_mx_c(x, n_bins, win_len, stride))

histMat_ref = zeros(n_bins, floor((length(x) - win_len + 1) / stride));

% Exhaustive computation for each window
tic
% Compute the edges for the histogram bins. Bin edges are set to be between
% 0 and 1 exactly here, 0 and 1 are included in the edges
histcounts_edges = (0:n_bins) / n_bins;
for i = 1:length(windows_loci_bt)
    idx = windows_loci_bt(i):(windows_loci_bt(i)+win_len-1);
    histMat_ref(:, i) = histcounts(x(idx), histcounts_edges);
end
toc

assert(isequal(histMat_bt, histMat_ref), 'Sliding window histograms do not match exhaustive computation.');
assert(isequal(histMat_custml, histMat_ref), 'Custom ML sliding window histograms do not match exhaustive computation.');
assert(isequal(histMat_custmx, histMat_ref), 'Custom MX sliding window histograms do not match exhaustive computation.');
assert(isequal(histMat_fullmx, histMat_ref), 'Full MX sliding window histograms do not match exhaustive computation.');
assert(isequal(histMat_mx_c, histMat_ref), 'MEX C sliding window histograms do not match exhaustive computation.');

assert(isequal(windows_loci_bt, windows_loci_custml), 'Window loci do not match between built-in and custom ML.');
assert(isequal(windows_loci_bt, windows_loci_custmx), 'Window loci do not match between built-in and custom MX.');
assert(isequal(windows_loci_bt, windows_loci_fullmx), 'Window loci do not match between built-in and full MX.');
assert(isequal(windows_loci_bt, windows_loci_mx_c), 'Window loci do not match between built-in and MEX C.');

assert(isequal(edges_bt,histcounts_edges), 'Edges do not match between built-in and exhaustive computation.');
assert(isequal(edges_custml, histcounts_edges), 'Edges do not match between custom ML and exhaustive computation.');
assert(isequal(edges_custmx, histcounts_edges), 'Edges do not match between custom MX and exhaustive computation.');
assert(isequal(edges_fullmx, histcounts_edges), 'Edges do not match between full MX and exhaustive computation.');
assert(isequal(edges_mx_c, histcounts_edges), 'Edges do not match between MEX C and exhaustive computation.');

disp(['All tests in '  mfilename() ' passed successfully!']);
%------------- END OF CODE --------------

