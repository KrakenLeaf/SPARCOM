function X = FISTA(Y,A,B,M, N, Params)
%% Initialization
% Paramters
beta       = Params.Beta;
L          = Params.L0; 
lambda_bar = Params.LambdaBar;
lambda     = Params.Lambda;

% Init
t_k        = 1;
t_k_m1     = 1;
X          = zeros(M, N);
X_k_m1     = X;

%% Iterations
for (ii = 1:Params.IterMax)
    Z=X+((t_k_m1-1)/t_k)*(X-X_k_m1);
    grad_f_z=(A'*(A*(Z*B)-Y)*B');
    U=Z-(1/L)*grad_f_z;
    
    %disp(['Iteration #:' num2str(ii)]);
    
    X_k_m1=X;
    X=sign(U).*(abs(U)-lambda/L).*(abs(U)>lambda/L);
    
    % Projection onto the non-negative orthant
    if (Params.NonNegOrth)
        X = real(X);                         % Orthogonal projection from {\mathbb C}^n to {\mathbb R}^n
        X = max(X,0);                        % Orthogonal projection from {\mathbb R}^n to {\mathbb R}^n_{+}
    end
    
    % Step update
    t_k_m1=t_k;
    t_k=(1+sqrt(4*t_k*t_k+1))/2;
    lambda=max(beta*lambda,lambda_bar);
end


