function [A_hat, E_hat, iter] = rpca_tv(M, lambda, tol, maxIter)

% Oct 2009
% This matlab code implements the inexact augmented Lagrange multiplier
% method for Robust PCA.
%
% D - m x n matrix of observations/data (required input)
%
% lambda - weight on sparse error term in the cost function
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
%
% maxIter - maximum number of iterations
%         - DEFAULT 1000, if omitted or -1.
%
% Initialize A,E,Y,u
% while ~converged
%   minimize (inexactly, update A and E only once)
%     L(A,E,Y,u) = |A|_* + lambda * |E|_1 + <Y,D-A-E> + mu/2 * |D-A-E|_F^2;
%   Y = Y + \mu * (D - A - E);
%   \mu = \rho * \mu;
% end
%
% Minming Chen, October 2009. Questions? v-minmch@microsoft.com ;
% Arvind Ganesh (abalasu2@illinois.edu)
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing

% addpath PROPACK;
%
[m, n] = size(M);

if nargin < 2
    lambda = 1 / sqrt(max(m,n));
end

if nargin < 3
    tol = 1e-6;
elseif tol == -1
    tol = 1e-6;
end

if nargin < 4
    maxIter = 100;
elseif maxIter == -1
    maxIter = 100;
end

% initialize
Y = M;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

A_hat = zeros( m, n);
% A_hat=M;
E_hat = zeros( m, n);
mu = 1.25/norm_two % this one can be tuned
mu_bar = mu * 1e7
rho = 1.5          % this one can be tuned
d_norm = norm(M, 'fro');


iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = 5;
while ~converged
    iter = iter + 1;
    
        temp_T = M - A_hat + (1/mu)*Y;
    %     E_hat = max(temp_T - lambda/mu, 0);
    %     E_hat = E_hat+min(temp_T + lambda/mu, 0);
    E_hat=proxTV(temp_T,lambda/mu);
    
    if choosvd(n, sv) == 1
        [U S V] = lansvd(M - E_hat + (1/mu)*Y, sv, 'L');
    else
        [U S V] = svd(M - E_hat + (1/mu)*Y, 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/mu));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
    
    A_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';

    
    total_svd = total_svd + 1;
    
    Z = M - A_hat - E_hat;
    
    Y = Y + mu*Z;
    mu = min(mu*rho, mu_bar);
    
    %% stop Criterion
    stopCriterion = norm(Z, 'fro') / d_norm;
    if stopCriterion < tol
        converged = true;
    end
    
    if mod( total_svd, 5) == 0
        disp(['#svd ' num2str(total_svd) ' r(A) ' num2str(rank(A_hat))...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);
    end
    
    if ~converged && iter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;
    end
end
end
    function u=proxTV(g,lbd)
    tau=0.249;

    udual=D(g);
    for i=1:16%8~16 are reasonable
        DDtz=D(Dt(udual)+g./lbd);
        udual=(udual-tau*DDtz)./(1+tau*abs(DDtz));
    end
    u=g + Dt(udual.*lbd);  
    end
function res=D(x)
res=x(:,[2:end,end])-x;
end

function res=Dt(x)
res=x(:,[1,1:end-1]) - x;
res(:,1)=-x(:,1);
res(:,end)=-x(:,end-1);
end