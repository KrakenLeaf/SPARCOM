function [ Xf ] = pfft2( X, M, iFlag, varargin )
% PFFT2 - Efficient partial FFT2 or IFFT2 on matix X. It holds that if F is an N X N Fourier matrix, then 
%
%         reshape(kron(F, F)*X(:), N, N) = F*X*F.' = fft2(X)
%
% We always assume that in the partial case, F is an M X N matrix with M < N. 
%  
% Syntax:
% -------
% [ Xf ] = pfft2( X, M, iFlag, varargin )
%
% Inputs:
% -------
% X        - Input 2D signal of arbitrary dimensions
% M        - Several cases:
%            iFLag == 'fft'  : M <= N (= length(X)), how many rows to take form a full N X N Fourier matrix. Relevant only to iFlag = 'fft'
%            iFlag == 'fft_t': M = N, "completion" to FFT of size N X N
%            iFlag == 'fft_h': M = N, "completion" to FFT of size N X N
% iFlag    - Type of calculation:
%            'fft'     : kron(F, F)*X(:) 
%            'fft_t'   : kron(F, F).'*X(:) (transpose)
%            'fft_h'   : kron(F, F)'*X(:)  (hermitian)
%            'fft_diag': Same as 'fft', only we assume that X is a zeros matrix except for one element in column specified by varargin{1}
% varargin - Only relevant for iFlag == 'fft_diag'
%            varargin{1}: Column number fot the non-zero entry in X
%            varargin{2}: We assume that X is an N X N matrix, where N = varargin{2}
%
% Outputs:  
% --------
% Xf    - Output signal
%
% Ver 1. Written by Oren Solomon, Technion I.I.T. 03-12-2015
%

switch lower(iFlag)
    case 'fft'
        %              Calculate: kron(F, F)*X(:)               %
        % ----------------------------------------------------- %
        if M > length(X); error('PFFT2: It must hold that M <= length(X)'); end
        
        % Temporary calculations
        X  = fft(X);
        X  = fft(X(1:M, :).');
        
        % Output
        Xf = X(1:M, :).';
    case 'fft_n_mid'
        MidRange = varargin{1};
        
        % Temporary calculations
        X  = fft(X);
        X  = fft([X(1:MidRange + mod(M,2),:); X(end - MidRange + 1:end,:)].');
        
        % Output
        Xf = [X(1:MidRange + mod(M,2),:); X(end - MidRange + 1:end,:)].';
    case 'fpft'
        %              Calculate: kron(F, F)*X(:)               %
        % ----------------------------------------------------- %
        % This is equivalent to 'fft', only with a different realization - Works
        % only when M and N are powers of 2!
        N = varargin{1};
        X = fpft(X, N, M);
        
        % Output
        Xf = fpft(X.', N, M).';
    case 'fft_diag'
        % We assume that X is a zeros matrix except for one element in column ii
        % Additional inputs
        ii = varargin{1};
        N = varargin{2};
        
        % fft only on relevant column
        Col_ii = fft(X(:, ii));

        % fft on padded matrix
        Xf = fft(padarray(Col_ii(1:M).', ii - 1, 0, 'pre'), N, 1);
        
        % Output
        Xf = Xf(1:M, :).';
    case 'fft_t'
        %              Calculate: kron(F, F).'*X(:)             %
        % ----------------------------------------------------- %       
        % Stage 1: Z = F.'*X      
        X = fft(X, M);
        
        % Stage 2: Xf = Z*F
        Xf = transpose(fft(X.', M));
    case 'fft_h'
        %              Calculate: kron(F, F)'*X(:)              %
        % ----------------------------------------------------- %   
        % Stage 1: Z = F'*X
        X = M*ifft(X, M);               % Factor M is needed for x = pfft{fft{x}} (M is N)
        
        % Stage 2: Xf = Z*conj(F)
        Xf = ctranspose(fft(X', M));    % Conjugate transpose on the fft
    case 'fft_h_n_mid'
        MidRange = varargin{1};         % M is actually N!!!
        
        k = (0:M - 1)';
        
        % Stage 1: Z = F'*X             %
        X = M*bsxfun(@times, ifft(ifftshift(X, 1), M), exp(1j*2*pi*k*MidRange/2/M));
        
        % Stage 2: Xf = Z*conj(F)
        Xf = ctranspose( bsxfun(@times, fft(fftshift(X', 1), M), exp(1j*2*pi*k*MidRange/2/M)) );
    otherwise
        error('PFFT2: iFlag not supported.');
end










