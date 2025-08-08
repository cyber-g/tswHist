function [histMat,strided_windows_loci,edges] = tswHist(input_norm, n_bins, win_len, stride, variant)
% TSWHIST - Fast sliding window histogram computation for 1D signals.
%   [histMat, strided_windows_loci, edges] = tswHist(input_norm, n_bins, win_len, stride, variant)
%
%   Computes histograms over sliding windows using efficient differential
%   updates. The local hist_int matlab function is adapted from the core of the
%   version found at: https://github.com/cyber-g/FastHist (only the essential
%   logic is retained). The hist_int_mx MEX function is an exact C transcription
%   of this local hist_int.
% 
%   Use tswHist_mx for a faster implementation of this function using MEX.
%
%   Inputs:
%     input_norm - Normalized input vector (1D signal) (in [0,1])
%     n_bins     - Number of histogram bins (integer > 2)
%     win_len    - Sliding window length
%     stride     - Stride for sliding window (default: 1)
%     variant    - 'builtin', 'custom-ml', or 'custom-mx' (default: 'builtin')
%
%   Outputs:
%     histMat              - n_bins x num_windows matrix of histograms
%     strided_windows_loci - Start indices of each window
%     edges                - Bin edges used for histogramming
%
%   See also: histcounts, hist_int_mx
%
%   Project: tswHist (https://github.com/cyber-g/tswHist)
%
%   Related paper: 
%       S. Perreault and P. Hebert, "Median Filtering in Constant Time," in IEEE
%       Transactions on Image Processing, vol. 16, no. 9, pp. 2389-2394, Sept.
%       2007, doi: 10.1109/TIP.2007.902329. URL:
%       https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4287006&isnumber=4286981
% 
%   Thanks:
%       Denis Gilbert (2025). M-file Header Template
%       (https://www.mathworks.com/matlabcentral/fileexchange/4908-m-file-header-template),
%       MATLAB Central File Exchange. Accessed on April 26, 2025.
%
%   License: GNU General Public License v3.0
%
%   Author: Germain PHAM
%   C2S, Télécom Paris, IP Paris
%   August 2025; Last revision: 

%------------- BEGIN CODE --------------

    % Validate inputs
    if nargin < 3
        error('Not enough input arguments. Provide input vector, window length, and stride (optional).');
    end
    if nargin < 4
        stride = 1;
    end
    if nargin < 5
        variant = 'builtin';
    end
    if ~isvector(input_norm)
        error('Input must be a vector.');
    end
    % n_bins must be a positive integer larger than 2
    if ~isscalar(n_bins) || n_bins <= 2 || n_bins ~= floor(n_bins)
        error('Number of bins must be a positive integer scalar larger than 2.');
    end
    if ~isscalar(win_len) || ~isscalar(stride)
        error('Window length and stride must be scalars.');
    end
    % Check that stride is smaller than win_len to guarantee overlap ; No
    % overlap makes no sense here because histogram are computed using a
    % differential approach
    if stride >= win_len
        error('Stride must be less than window length to ensure overlap.');
    end

    % Compute the strided windows loci
    strided_windows_loci = 1:stride:(length(input_norm) - win_len + 1);

    % Compute the number of strided windows
    num_strided_windows  = floor(length(strided_windows_loci));

    % Initialize the histogram matrix
    histMat              = zeros(n_bins, num_strided_windows);

    % The normalization is left outside this function for more flexibility The
    % input vector is expected to be included in [0,1] (not necessarily exactly
    % occupying this range)
    input_int  = floor(input_norm * n_bins); 

    % Compute the histogram for the first window
    switch variant
        case 'builtin'
            % Use the built-in histcounts function
            bufferHist = histcounts(input_int(1:win_len), 0:n_bins);
        case 'custom-ml'
            % Use a personal open source implementation of histogram
            % Patch the input vector to match histcounts behavior
            input_int(input_int==n_bins) = n_bins-1; 
            bufferHist = hist_int(input_int(1:win_len), n_bins);
        case 'custom-mx'
            % Use a personal open source implementation of histogram
            % Patch the input vector to match histcounts behavior
            input_int(input_int==n_bins) = n_bins-1; 
            bufferHist = hist_int_mx(input_int(1:win_len), n_bins);
        case 'pushHist'
            % Patch the input vector to match histcounts behavior
            input_int(input_int==n_bins) = n_bins-1; 
            bufferHist = pushHist(zeros(1,n_bins),input_int(1:win_len));
        otherwise
            error('Unknown variant specified. Use "builtin", "custom-ml", or "custom-mx".');
    end
    histMat(:, 1) = bufferHist;

    % Patch the input vector so as to match histcounts behavior in the case of hist_int
    % https://www.mathworks.com/help/matlab/ref/double.histcounts.html
    % Each bin includes the local leading edge, but does not include the local
    % trailing edge, except for the last bin which includes both edges.
    input_int(input_int==n_bins) = n_bins-1; % Ensure the max value is integrated to the last bin

    

    % Offsets for the popping and pushing elements indices
    io_offsets = (-(stride-1):0); 

    % now for each subsequent window, 
    for i = 2:length(strided_windows_loci)
        pop_out_indexs = strided_windows_loci(i)   - 1      + io_offsets;
        push_in_indexs = strided_windows_loci(i)+(win_len-1)+ io_offsets;
        
        % Compute differential histograms
        bufferHist  = popHist(bufferHist, input_int(pop_out_indexs));
        bufferHist  = pushHist(bufferHist, input_int(push_in_indexs));

        % Store the histogram
        histMat(:, i) = bufferHist;

    end
    % Compute the edges for the histogram bins. Bin edges are set to be between
    % 0 and 1 exactly here, 0 and 1 are included in the edges
    edges = (0:n_bins)/n_bins;

end

function bufferHist = pushHist(bufferHist, input_int)
    % Increment the histogram counts for the new elements
    for j = 1:length(input_int)
        % Ignore bounding test
        bufferHist(input_int(j)+1) = bufferHist(input_int(j)+1) + 1;
    end
end

function bufferHist = popHist(bufferHist, input_int)
    % Decrement the histogram counts for the removed elements
    for j = 1:length(input_int)
        % Ignore bounding test
        bufferHist(input_int(j)+1) = bufferHist(input_int(j)+1) - 1;
    end
end

function hist_vec = hist_int(input_int, n_bins)
    % Custom histogram implementation for integer input
    % Please note that computing the histogram for a vector of integers from
    % scratch is equal to use pushHist on an empty histogram vector with a full
    % input vector
    hist_vec = zeros(1, n_bins);
    hist_vec = pushHist(hist_vec, input_int);
end

%------------- END OF CODE --------------