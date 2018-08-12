function [ X_out, S_out ] = MatOMP( Y, A, B, Params )
%MATOMP performs the OMP algorithm in matrix form
%
% Syntax:
% -------
% [ X_out, S_out ] = MatOMP( Y, A, B, Params )
%
% Inputs:
% -------
% Y      - Inputs data matrix
% A      - Input sensing matrix A
% B      - Input sensing matrix B
% Params - Algorithm parameters
%          'Tol'        - Stopping criterion tolerance
%          'IterMax'    - Maximum number of iterations
%          'NonNegOrth' - Projection onto the non-negative orthant flag (0 - no, 1 - yes)
%
% Outputs:
% --------
% X_out  - Output solution
% S_out  - Output support
%

%% Initialization
method = 1;
C = [];

R = Y;                                                  % Residual
S = [];                                                 % empty index set

kk = 1;                                                 % Iteration index

[MA, NA] = size(A);
[MB, NB] = size(B);

X = zeros(NA, NB);

% Exit the algorithm if the specified tolerance is negative
if (Params.Tol < 0)
    X_out = X;
    S_out = [];
    
    return;
end

% Normalize matrices A and B
[ A, Da, Dainv ] = NormalizeMat( A );
[ B, Db, Dbinv ] = NormalizeMat( B );

frob =  norm(Y,'fro');                                  % Initialization of the tolerance parameter
frob_prev = 0;
difsol = abs(frob - frob_prev);

%% Iterations
while (kk <= Params.IterMax & difsol >  Params.Tol)
    if (mod(kk,10) == 0)
        disp(['Iteration #' num2str(kk) ' : Frobenius norm diff is ' num2str(difsol)]);
    end
    
%     clear x;
    frob_prev = frob;
    
    % Sweep stage
    SweepStack = abs(B'*R'*A);
    [~, Icol] = max(max(SweepStack));                   % Find column of max value
    [~, Irow] = max(SweepStack(:,Icol));                % Find row of max value
    
    % Update 2D support
    S = [S; [Icol, Irow]];
    
    % Breaking condition
    if (kk > 1)
        Iuni1 = find(S(1:end-1,1) == Icol);
        if (~isempty(Iuni1))
            Iuni2 = find(S(Iuni1,2) == Irow);
            if (~isempty(Iuni2))
                break;
            end
        end
    end
    
%     % Find if the current support estimation waas already found or not - it might happen that the same indices are rechosen.
%     % In this case, eliminate that entry in SweepStack.
%     Iuni1 = find(S(:,1) == Icol);
%     Iuni2 = find(S(:,2) == Irow);
%     while (length(Iuni1) > 1 & length(Iuni2) > 1)
%         % Eliminate recurrent element in SweepStack
%         SweepStack(Irow, Icol) = 0;
%     
%         % Perform sweep again
%         [~, Icol] = max(max(SweepStack));                   % Find column of max value
%         [~, Irow] = max(SweepStack(:,Icol));                % Find row of max value
%     
%         % Reupdate 2D support
%         S = [S(1:end-1,:); [Icol, Irow]];
%     
%         Iuni1 = find(S(:,1) == Icol);
%         Iuni2 = find(S(:,2) == Irow);
%     end

    % Choose solution method - actually they produce the same result.
    if (method == 1)
        % New signal estimate
        d_temp = B(:,S(:,2))'*Y';                           % Create the vector d
        d = zeros(kk,1);
        for (ii = 1:kk)
            d(ii) = d_temp(ii,:)*A(:,S(ii,1));
        end
        
        D = zeros(kk,kk);                                   % Create the matrix D
        for (ii = 1:kk)
            for (jj = 1:kk)
                D(ii,jj) = B(:,S(jj,2))'*B(:,S(ii,2))*A(:,S(ii,1))'*A(:,S(jj,1));
            end
        end
        
        %     x = pinv(D)*d;
        x = pinv(conj(D))*d;                                % Vector estimate, whose size increases by 1 in each iteration
        %     x = D\d;
        
        % Projection onto the non-negative orthant
        if (Params.NonNegOrth)
            x = real(x);                                    % Orthogonal projectino from {\mathbb C}^n to {\mathbb R}^n
            x = max(x,0);                                   % Orthogonal projection from {\mathbb R}^n to {\mathbb R}^n_{+}
        end
        
        % New approximation
        Q = zeros(MA, MB);
        for (ii = 1:kk)
            Q = Q + x(ii)*A(:,S(ii,1))*B(:,S(ii,2))';
        end
        
        % Compute the residual
        R = Y - Q;
    else
        C = [C vec(A(:,S(end,1))*B(:,S(end,2))')];
        x = pinv(C)*Y(:);
        
        % Projection onto the non-negative orthant
        if (Params.NonNegOrth)
            x = real(x);                                    % Orthogonal projectino from {\mathbb C}^n to {\mathbb R}^n
            x = max(x,0);                                   % Orthogonal projection from {\mathbb R}^n to {\mathbb R}^n_{+}
        end
        
        R = reshape(Y(:) - C*x, MA, MB);
    end

    
    % Advance iteration
    kk = kk + 1;
    
    % Frobenius norm stopping criterion
    for (ii = 1:kk-1)
        X(S(ii,1),S(ii,2)) = x(ii);
    end
    frob   = norm(Y - A*(X*B'),'fro');                    % Frobenius norm
    difsol = abs(frob - frob_prev);
end

%% Output
for (ii = 1:kk-1)
    X(S(ii,1),S(ii,2)) = x(ii);
end

% X_out = X;
X_out = Da*(X*Db);
S_out = S(1:kk-1,:);

% disp(['Number of iterations performed by MatOMP: ' num2str(kk-1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          Ancillary functions                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Xnrm, D, Dinv ] = NormalizeMat( X )
% L2 Normalilze the columns of matrix X

[~, N] = size(X);
D = zeros(N, N);
Dinv = zeros(N, N);
for (ii = 1:N)
	D(ii, ii) = 1/norm(X(:,ii));
    Dinv(ii, ii) = norm(X(:,ii));
end
Xnrm = X*D;

%% Obsolete Code Segment - instead of "Breaking condition"

%     % Breaking condition
%     Iuni1 = find(S(:,1) == Icol);
%     Iuni2 = find(S(:,2) == Irow);
%     if (length(Iuni1) > 1 && length(Iuni2) > 1) break; end;


%         % Find if the current support estimation waas already found or not - it might happen that the same indices are rechosen.
%         % In this case, eliminate that entry in SweepStack.
%         Iuni1 = find(S(:,1) == Icol);
%         Iuni2 = find(S(:,2) == Irow);
%         while (length(Iuni1) > 1 & length(Iuni2) > 1)
%             % Eliminate recurrent element in SweepStack
%             SweepStack(Irow, Icol) = 0;
%     
%             % Perform sweep again
%             [~, Icol] = max(max(SweepStack));                   % Find column of max value
%             [~, Irow] = max(SweepStack(:,Icol));                % Find row of max value
%     
%             % Reupdate 2D support
%             S = [S(1:end-1,:); [Icol, Irow]];
%     
%             Iuni1 = find(S(:,1) == Icol);
%             Iuni2 = find(S(:,2) == Irow);
%         end






















