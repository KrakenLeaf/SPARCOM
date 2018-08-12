function [ x, X_stack ] = pFISTA_Iterative( Y, H, Params )
%FPG_ITERATIVE - Iterative weighted l1 minimization, combined with a FISTA type algorithm
% 
%              min_{x} \lambda|| diag{e}x ||_1 + 0.5|| Y - AXA^H ||_2^2
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
% [ x, X_stack ] = pFISTA_Iterative( Y, H, Params )
%
% Inputs:
% -------
% Y                 - [My, Ny, Ky] cube of 2D frames
% A                 - [Ma, Na] Left measurement matrix   
% B                 - [Mb, Nb] Right measurement matrix
% W                 - Weights vector of lemgth Ky
% Params.L          - Lipschitz constant of f(X) = \sum_{j=1}^{m}w_j||Y_j -
%                     \sum_{i=1}^{n}M_i*x_{ij}||_F^2. If empty, the function calculates it herself
% Params.Lambda     - Regularization parameter. Larger values lead to sparser solutions
% Params.IterMax    - Maximum number of iterations 
% Params.LoadConsts - 0 - calculate Mtr and V, 1 - load them from memory (since for large arrays, the calculation is time consuming)
% Params.ConstName  - Name of file to load, which contains Mtr and V  
% Params.SaveConsts - 0 - do not save Mtr and V, 1 - save them
% Params.NumOfSteps - NUmber of outer-loop iterations
% Params.eps        - Epsilon for the l1 weights
%
% Outputs:
% --------
% X                 - Solution
% L                 - Lipschitz constant used in the calculation
%
% Written by Oren Solomon, Technion I.I.T. Ver 1. 28-12-2015
%

global VERBOSE

%% Initialization
% -----------------------------------------------------------------------------------------------------------------------
TotIter = Params.NumOfSteps;                            % Total number of "enveloping" iterations
eps     = Params.eps;                                   % Regularization parameter of the "enveloping" iterative scheme
N       = Params.N;                                     % Size of x
L       = Params.L0;                                    % Lipschitz constant

% Initial weights
e = ones(N, 1);

X_stack = [];

try
    LargeScaleFlag = Params.LargeScale;
catch
    LargeScaleFlag = 0;              % Default value - small scale
end

%% Calculate constants
% -----------------------------------------------------------------------------------------------------------------------
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

%% Iterations
% -----------------------------------------------------------------------------------------------------------------------
for ii = 1:TotIter 
    % Output to screen
    if VERBOSE; disp(['pFISTA_Iterative scheme: Iteration #' num2str(ii)]); disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'); end;
    
    % Perform FPG algorithm
    [ x ] = pFISTA_diag_it( e, v, S, L, Params );
    
    % Update weights      
    e = 1./(abs(x) + eps);                      % E is a vector - This is an l1 norm for each row of X
    
    % Stack solutions
    X_stack = [X_stack x];
end
if VERBOSE; disp('Done performing pFISTA_Iterative.'); end;

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

%% An auxiliary function for v_calc, to speed up computations - output is a scalar - large data version
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
Q = Nsqrt*ifft(t, Nsqrt); 

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
Q = Nsqrt*ifft(t, Nsqrt); 

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
S = fft2(reshape(A1, sqrt(N), sqrt(N)));

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

X = arrayfun(@(ii) FirstElement(pfft2(reshape(full(Y(:, ii)), sqrt(My), sqrt(My)), sqrt(N), 'fft_h')), 1:Ny) ;

%% Take first element of a matrix
function a = FirstElement(Q)
a = Q(1, 1);

%% Vectorize a matrix
function v = vec( x )
v = x(:);


