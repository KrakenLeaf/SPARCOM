function [ MatCorr ] = CalcCorr( MatIn, Dims, parFlag )
%CALCCORR - Calculate correlation of the input matrix
%   We assume the following structure on MatIn: Row index indicate element in some
%   vector y(t), column index indicate time-stamp.
%
% Syntax:
% -------
% [ MatCorr ] = CalcCorr( MatIn, varargin )
%
% Inputs:
% -------
% MatIn       - Input matrix. Each column is a vector of pixels in a
%               different time-stamp (Pixel number X Time-stamp: M X N)
% varargin{1} - Row dimension (Optional)
% varargin{2} - Column dimension (Optional)
% varargin{3} - Perform parallel computation (Optional)
%
% Outputs:
% --------
% MatCorr     - A matrix of the size M^2 X N. Each column is a vectorization of the
%               correlation matrix of the columns of MatIn, for a different time-lag (Tau = 0,1,...,N-1).
%
% Ver 1: Written by Oren Solomon, 03-08-2015
%

% Determine MatIn dimensions
if ~isempty(Dims)
    M = Dims(1);
    N = Dims(2);
else
    [M, N] = size(MatIn);
end

% Initialize
MatCorr = zeros(M^2, N);

% Calculate mean - for each row (each pixel's time-trace)
CalcBias = mean(MatIn, 2);

% Perform correlation calculations
% ----------------------------------------------------------------------------------------------------------------------------------------------
if parFlag
    % Perform parallel computation
    % ----------------------------
    parfor Increment = 0:(N-1)
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
%             TmpMat = TmpMat + MatIn(:, ii)*(MatIn(:, jj))';     % No mean subtraction
%             TmpMat = TmpMat + (MatIn(:, ii) - mean(MatIn(:, ii)))*(MatIn(:, jj) - mean(MatIn(:, jj)))';     % With mean subtraction
            TmpMat = TmpMat + (MatIn(:, ii) - CalcBias)*(MatIn(:, jj) - CalcBias)';     % With mean subtraction
        end
        
        % Store calculation in vector form - Each column is an average
        MatCorr(:, Increment + 1) = (1/(N - Increment))*reshape(TmpMat, M^2, 1);
    end
else
    % Perform serial computation
    % --------------------------
    for Increment = 0:(N-1)
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
%             TmpMat = TmpMat + MatIn(:, ii)*(MatIn(:, jj))';                             % No mean subtraction
            TmpMat = TmpMat + (MatIn(:, ii) - CalcBias)*(MatIn(:, jj) - CalcBias)';     % With mean subtraction
        end
        
        % Store calculation in vector form - Each column is an average
        MatCorr(:, Increment + 1) = (1/(N - Increment))*reshape(TmpMat, M^2, 1);
    end
end

% Take roughly 90% of the correlations
% MatCorr = MatCorr(:,1:floor(0.9*N));







