function [ x, lambda_stack ] = pFISTA_diag( Y, H, Params )
% PFISTA_DIAG - Performs the Fast Proximal Gradient method on the problem
% 
%              min_{x} \lambda|| x ||_1 + 0.5|| Y - AXA^H ||_F^2
%   
% with X being diagonal with a sparse diagonal x. A is assumed to have the following 
% specific structure: A = H*kron(Fp, Fp), with H being a diagonal matrix  and Fp 
% is a partial Fourier matrix. This structure leads to a very efficient implementation 
% (in terms of run-time and memory usage) is using fft and ifft operations. 
%
% It is also assumed that x is non-negative and real (which is true for a correlation matrix, for instance)
%
% Syntax:
% -------
% [ x ] = pFISTA_diag( Y, H, Params )
%
% Inputs:
% -------
% Y      - Input M^2 X M^2 observations matrix (Generaly complex)
% H      - M^2 X M^2 diagonal matrix (should be in sparse format)
% Params - Additional algorithmic parameters
%          Beta      : Regularization parameter (0,1) - best to take very close to 1
%          L0        : Lipschitz constant (norm(kron(A, A), 2))
%          LambdaBar : Regularization parameter. Should be very small, e.g. 1e-7
%          Lambda    : Sparsity regularization parameter. Value depends on the problem
%          N         : Reconstructed Xr is of size N^2 X N^2
%          isPos     : 1 if X entries are supposed to be non negative
%          IterMax   : Maximum number of iterations
%          LargeScale: 1 if N > 500^2 (roughly). Computation slows dramatically, but can handle very large datasets
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
L          = Params.L0;              % Lipschitz constant
lambda_bar = Params.LambdaBar;
lambda     = Params.Lambda;          % l1 regularization parameters
N          = Params.N;               % Length of x (Actually N^2)

M2         = size(Y, 1);

% Init
t          = 1;
t_prev     = 1;

% Memory allocation
x          = zeros(N, 1);
x_prev     = x;

try
    LargeScaleFlag = Params.LargeScale;
catch
    LargeScaleFlag = 0;              % Default value - small scale
end

if VERBOSE; disp('Performing constant calculations...'); ConstTime = tic; end;
% Caclculate the vector v
if VERBOSE; fprintf('v_calc: ');t1 = tic; end;
if LargeScaleFlag == 1
    v = real(v_calc_large(Y, H, N)); % Real is to clean off any imaginary residuals. v should be real
else
    v = real(v_calc(Y, H, N));       % Real is to clean off any imaginary residuals. v should be real
end
if VERBOSE; toc(t1); end;

% Calculate S
if VERBOSE; fprintf('S_calc: ');t2 = tic; end;
S = S_calc(H, N);
if VERBOSE; toc(t2); end;

if isempty(L)
    L = max(max(S));
end

if VERBOSE; disp('done.'); toc(ConstTime); end;

% Automatic regularization parameter update
if isempty(lambda)
    lambda = 0;
    AutoLambdaFlag = 1;
else
    AutoLambdaFlag = 0;
end

% lambda_stack = zeros(1, Params.IterMax);
lambda_stack = [];

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
    x = sign(g).*max(abs(g) - lambda/L, 0);
    
    % Projection onto the non-negative orthant, only for a non-negative constraint
    if Params.NonNegOrth == 1
        x(x < 0) = 0;
    end
    
    % Fast exit
    if norm(x) == 0
        return;
    end
    
    % Parameter updates for next iteration
    t_prev = t;
    t = 0.5*(1 + sqrt(4*t^2 + 1));
    
    % Automatic lambda update
    if AutoLambdaFlag
        Tmp     = LAH(ctranspose(LAH(diag(x), H)), H);
        lambda  = (norm(Y - ctranspose(Tmp), 'fro')^2)/(N^2 - M2); %(M2 - N^2) 
    else
        lambda = max(beta*lambda, lambda_bar);
    end
    
    %lambda_stack(kk) = lambda;
    % Speed-up: Check if the lambda difference is smaller than a predefined threshold
    lambda_stack = [lambda_stack lambda];
%     if kk > 2 && AutoLambdaFlag && abs(lambda_stack(end) - lambda_stack(end - 1)) < 1e-13 %1e-12 % This threshold parameter can be tuned
%         return 
%     end
    
    if VERBOSE; toc(Titer); end;
end
if VERBOSE; disp('Done pFISTA_diag.'); end;

%% Auxiliary functions
% -----------------------------------------------------------------------------------------------------------------------
%% Calculate the constant v - large data version
function v = v_calc_large(Y, H, N)
% Init
[~, M] = size(H);
Nsqrt = sqrt(N); 
Msqrt = sqrt(M);

P = H' * (Y * H);     % M^2 X M^2

v = zeros(N, 1);
for ii = 1:N
    qii = v_calc_inner_loop_large(ii, P, N, M);    % M X 1 vector
    v(ii) = lkfft(ii, qii, Nsqrt, Msqrt);
end
% v = arrayfun(@(ii) conj( lkfft(ii, v_calc_inner_loop(ii, P, N, M), Nsqrt, Msqrt ) ), (1:N)'); % Same runtime as using the FOR loop

%% An auxiliary function for v_calc, to speed up computations - output is a scalar - large data version. Caution: very slow
function vt = v_calc_inner_loop_large(ii, P, N, M)
Nsqrt = sqrt(N); 
Msqrt = sqrt(M);

vt = arrayfun(@(jj) lkfft(ii, P(:, jj), Nsqrt, Msqrt), (1:M)');     % M X 1 vector

%% Calculate the constant v
function v = v_calc(Y, H, N)
% Init
[M, M] = size(H);
v = zeros(N, 1);
Nsqrt = sqrt(N);
Msqrt = sqrt(M);

P = H'*LAH_H(Y, H, N)';            % M X N: P = H^H*Z, Z^H = A^H*Y

v = arrayfun(@(ii) v_calc_inner_loop(ii, P(:, ii), Nsqrt, Msqrt), (1:N)' );

%% An auxiliary function for v_calc, to speed up computations - output is a scalar
function vt = v_calc_inner_loop(ii, Pv, Nsqrt, Msqrt)
% Step 0: Determine indices
ki = floor((ii - 1)/Nsqrt) + 1;
li = mod(ii, Nsqrt) + Nsqrt*(mod(ii, Nsqrt) == 0);

% Step 1: Convert to M X M matrix
t = reshape(Pv, Msqrt, Msqrt);

% Step 2: Q is an N X M matrix
Q = Nsqrt*ifft(t, Nsqrt); % N*ifft(t, Nsqrt); Why Nsqrt^2 and not Nsqrt ?

% Step 3: Output is an N X 1 vector
q2 = fft(Q(li, :)', Nsqrt);

% Output
vt = q2(ki);

%% Implementation of a_i^H*p
function a = lkfft(ii, Pv, Nsqrt, Msqrt)
% Step 0: Determine indices
ki = floor((ii - 1)/Nsqrt) + 1;
li = mod(ii, Nsqrt) + Nsqrt*(mod(ii, Nsqrt) == 0);

% Step 1: Convert to M X M matrix
t = reshape(Pv, Msqrt, Msqrt);

% Step 2: Q is an N X M matrix
Q = Nsqrt*ifft(t, Nsqrt); % N*ifft(t, Nsqrt); Why Nsqrt^2 and not Nsqrt ?

% Step 3: Output is an N X 1 vector
q2 = fft(Q(li, :)', Nsqrt);

% Output
a = q2(ki);

%% Calculate A1 which is needed for the gradient calculation 
function S = S_calc(H, N)
% Step 0
H2  = diag( abs(diag(H)).^2 );

% Step 1: M^2 X 1 vector q
q  = ctranspose( LAH_I(H2, N) );

% Step 2: N^2 X 1 vector A1
A1 = abs(LAH_H(q, 1, N)).^2;

% Step 3: Calculate eigenvalues sqrt(N) X sqrt(N) matrix (N eigenvalues)
S = fft2(reshape(A1, sqrt(N), sqrt(N)));        % 1/N ? - does not seem to affect reconstruction performance


function g = AHA(x, S, N)
g = cell2mat( arrayfun(@(ii)  ifft2(S .* fft2( reshape(x(:, ii), sqrt(N), sqrt(N)) )), 1:N, 'UniformOutput', false));
%B = ifft2(S .* fft2( reshape(x, sqrt(N), sqrt(N)) ));
%g = real(B(:));

%% Calculate the gradient step: \nabla f(x) = |A^H*A|^2*x - v. |A^H*A|^2 is a BCCB matrix and admits a fast matrix-vector multiplication
function g = Grad(x, S, v, N)
B = ifft2(S .* fft2( reshape(x, sqrt(N), sqrt(N)) ));
g = real(B(:)) - v;

%% Left AH: Implementation of A*Y = H*kron(F, F)*Y efficiently, using FFT operations
function X = LAH(Y, H)
% Determine dimensions
[My, Ny] = size(Y);
[Mh, Nh] = size(H);

% X = H * vec(pfft2(reshape(Y, sqrt(My), sqrt(My)), sqrt(Mh), 'fft'));
X = H * cell2mat( arrayfun(@(ii) vec(pfft2(reshape(full(Y(:, ii)), sqrt(My), sqrt(My)), sqrt(Mh), 'fft')), 1:Ny, 'UniformOutput', false) );

%% Left AH Hermitian: Implementation of A^H*Y = kron(F, F)^H*H^H*Y efficiently, using FFT operations
function X = LAH_H(Y, H, N)
% Determine dimensions
[My, Ny] = size(Y);

Z = H' * Y;

X = cell2mat( arrayfun(@(ii) vec(pfft2(reshape(Z(:, ii), sqrt(My), sqrt(My)), sqrt(N), 'fft_h')), 1:Ny, 'UniformOutput', false) );

%% Similar to LAH_H, only the output is a vector
function X = LAH_I(Y, N)
% Determine dimensions
[My, Ny] = size(Y);

X = arrayfun(@(ii) FirstElement(pfft2(reshape(full(Y(:, ii)), sqrt(My), sqrt(My)), sqrt(N), 'fft_h')), 1:Ny);

%% Take first element of a matrix
function a = FirstElement(Q)
a = Q(1, 1);

%% Vectorize a matrix
function v = vec( x )
v = x(:);
