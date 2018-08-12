function [ MatCorr ] = CalcCorr_Mat_V1( MatIn, Dims, MaxTimeLag, parFlag )
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

% Initialize
MatCorr = zeros(M, M, MaxTimeLag);

% Calculate mean - for each row (each pixel's time-trace)
CalcBias = mean(MatIn, 2);

% Subtruct the mean from the measurements
MatIn = MatIn - repmat(CalcBias, 1, N);

%% Perform correlation calculations
% ----------------------------------------------------------------------------------------------------------------------------------------------
if parFlag
    % Perform parallel computation
    % ----------------------------
    parfor Increment = 0:(MaxTimeLag - 1)
        % Initialize accumulator
        TmpMat = zeros(M, M);
        
        % --------------------------------------------------------------------------------------------------------------------------------------
        % Perform correlations - Note that there is redundancy for each time-lag:
        %           [ Tau = 0: E{y(t_1)y^H(t_1)}=...=E{y(t_N)y^H(t_N)}   ]
        %           [ Tau = 1: E{y(t_1)y^H(t_2)}=...=E{y(t_N-1)y^H(t_N)} ]
        %           [ ...                                                ]
        %           [ Tau = N: E{y(t_1)y^H(t_N)}                         ]
        % --------------------------------------------------------------------------------------------------------------------------------------
        for ii = 1:(N - Increment)
            jj = ii + Increment;
            TmpMat = TmpMat + MatIn(:, ii)*MatIn(:, jj)';     % With mean subtraction
        end
        
        % Store calculation in matrix form
        MatCorr(:, :, Increment + 1) = (1/(N - Increment))*TmpMat;
    end
else
    % Perform serial computation
    % --------------------------
    for Increment = 0:(MaxTimeLag - 1)
        % Initialize accumulator
        TmpMat = zeros(M, M);
        
        % --------------------------------------------------------------------------------------------------------------------------------------
        % Perform correlations - Note that there is redundancy for each time-lag:
        %           [ Tau = 0: E{y(t_1)y^H(t_1)}=...=E{y(t_N)y^H(t_N)}   ]
        %           [ Tau = 1: E{y(t_1)y^H(t_2)}=...=E{y(t_N-1)y^H(t_N)} ]
        %           [ ...                                                ]
        %           [ Tau = N: E{y(t_1)y^H(t_N)}                         ]
        % --------------------------------------------------------------------------------------------------------------------------------------
        for ii = 1:(N - Increment)
            jj = ii + Increment;
            TmpMat = TmpMat + MatIn(:, ii)*MatIn(:, jj)';     % With mean subtraction
        end
        
        % Store calculation in matrix form
        MatCorr(:, :, Increment + 1) = (1/(N - Increment))*TmpMat;
    end
end
