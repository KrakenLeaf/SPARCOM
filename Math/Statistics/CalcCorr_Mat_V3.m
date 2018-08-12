function [ MatCorr ] = CalcCorr_Mat_V3( MatIn, Dims, MaxTimeLag )
%CALCCORR - Calculate correlation of the input matrix
%   We assume the following structure on MatIn: Row index indicate element in some
%   vector y(t), column index indicate time-stamp.
%
% Syntax:
% -------
% [ MatCorr ] = CalcCorr_Mat( MatIn, Dims, NumTimeLags, parFlag )
%
% Inputs:
% -------
% MatIn      - Input matrix. Each column is a vector of pixels in a different time-stamp (Pixel number X Time-stamp: M X N)
% Dims       - [Row dimension, Column dimension] (Optional)
% MaxTimeLag - Maximum time-lag from zero (optional)
% parFlag    - Perform parallel computation (Optional)
%
% Outputs:
% --------
% MatCorr    - A cube of the size M X M X N. Each frame (3rd dimension) 
%              is the correlation matrix for a different time-lag (Tau = 0,1,...,N-1).
%
% Ver 1: Written by Oren Solomon, 03-08-2015 - Output is in vector form
% Ver 2: Written by Oren Solomon, 04-10-2015 - Output is in matrix form
% Ver 3: Written by Oren Solomon, 13-01-2016 - Efficient matrix calculation
%

%% Initialization
% ----------------------------------------------------------------------------------------------------------------------------------------------
% Determine MatIn dimensions
if ~isempty(Dims)
    M = Dims(1);
    N = Dims(2);
else
    [M, N] = size(MatIn);
end

if isempty(MaxTimeLag)
    MaxTimeLag = N;
end

% Calculate mean - for each row (each pixel's time-trace)
CalcBias = mean(MatIn, 2);

% Subtruct the mean from the measurements
MatIn = MatIn - repmat(CalcBias, 1, N);

%% Perform correlation calculations
% ----------------------------------------------------------------------------------------------------------------------------------------------
MatCorr = arrayfun(@(Increment) (1/(N - Increment))*(MatIn(:, 1:end - Increment)*MatIn(:, Increment + 1:end)'), 0:(MaxTimeLag - 1), 'UniformOutput', false);
MatCorr = cat(3, MatCorr{:});




































