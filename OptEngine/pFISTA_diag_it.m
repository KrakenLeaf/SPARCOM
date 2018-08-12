function [ x ] = pFISTA_diag_it( e, v, S, L, Params )
% PFISTA_DIAG - Performs the Fast Proximal Gradient method on the problem
% 
%              min_{x} \lambda|| diag{e}x ||_1 + 0.5|| Y - AXA^H ||_F^2
%   
% with X being diagonal with a sparse diagonal x. A is assumed to have the following 
% specific structure: A = H*kron(Fp, Fp), with H being a diagonal matrix  and Fp 
% is a partial Fourier matrix. This structure leads to a very efficient implementation 
% (in terms of run-time and memory usage) is using fft and ifft operations. 
%
% It is also assumed that x is non-negative and real (which is true for a correlation matrix, for instance)
%
% This pFISTA version has a diagonal weighting matrix diag{e} = E > 0 and is used inside pFISTA_Iterative.m.
%
% Syntax:
% -------
% [ x ] = pFISTA_diag_it( e, v, S, L, Params )
%
% Inputs:
% -------
% Y      - Input M^2 X M^2 observations matrix (Generaly complex)
% H      - M^2 X M^2 diagonal matrix (should be in sparse format)
% e      - N^2 X 1 weighting vector (non-negative entries)
% Params - Additional algorithmic parameters
%          Beta     : Regularization parameter (0,1) - best to take very close to 1
%          L0       : Lipschitz constant (norm(kron(A, A), 2))
%          LambdaBar: Regularization parameter. Should be very small, e.g. 1e-7
%          Lambda   : Sparsity regularization parameter. Value depends on the problem
%          N        : Reconstructed Xr is of size N^2 X N^2
%          isPos    : 1 if X entries are supposed to be non negative
%          IterMax  : Maximum number of iterations
%
% Output:
% -------
% x      - An estimate of the sparse N X 1 vector.
%
% Written by Oren Solomon, Technion I.I.T. Ver 1. 28-12-2015
%

global VERBOSE

%% Initializations
% -----------------------------------------------------------------------------------------------------------------------
% Paramters
beta       = Params.Beta;
lambda_bar = Params.LambdaBar;
lambda     = Params.Lambda;         % l1 regularization parameters
N          = Params.N;              % Length of x

% Init
t          = 1;
t_prev     = 1;

% Memory allocation
x          = zeros(N, 1);
x_prev     = x;

%% Iterations
% -----------------------------------------------------------------------------------------------------------------------
if VERBOSE; disp('pFISTA_diag: Running iterations...'); end;
for kk = 1:Params.IterMax
    Titer = tic;
    if VERBOSE; fprintf(['Iteration #' num2str(kk) ': ']); end;
    
    % Z update
    z = x + ((t_prev - 1)/t)*(x - x_prev);  
    
    % Gradient step: G = Z - (1/L)*( |A^H*A|^2*x - v )
    g = z - (1/L)*Grad(z, S, v, N);

    % Soft thresholding
    x_prev = x;
    x = sign(g).*max(abs(g) - e*lambda/L, 0);  % e is a vector 

    if Params.NonNegOrth == 1
       x(x < 0) = 0; 
    end
    
    % Parameter updates for next iteration
    t_prev = t;
    t = 0.5*(1 + sqrt(4*t^2 + 1));
    
    lambda = max(beta*lambda, lambda_bar);
    
    if VERBOSE; toc(Titer); end;
end
if VERBOSE; disp('done.'); end;

%% Auxiliary functions
% -----------------------------------------------------------------------------------------------------------------------
%% Calculate the gradient step: \nabla f(x) = |A^H*A|^2*x - v. |A^H*A|^2 is a BCCB matrix and admits a fast matrix-vector multiplication
function g = Grad(x, S, v, N)
B = ifft2(S .* fft2( reshape(x, sqrt(N), sqrt(N)) ));
g = real(B(:)) - v;
