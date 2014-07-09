function [ F, F_star, s] = RSGNM(W, K, lambda, beta, tau, tol, n_iter, F_init)
% RSRGM is the main algorithm described in Algorthm 1. Instead of initializing propensity matrix F based on 
% the resulting complexes of SPICi[16] as described in Algorithm 1, here the code initializes F as a matrix of 
% n by k by random. You can also initialize it as the way presented in Step
% 1),2) and 3), and  input it as the value of parameter 'F_init'.
%
% Inputs:
%   W: adjacent matrix of PPI network under consideration.
%   K: maximum number of possible protein complexes.
%   lambda: penalization parameter for smooth regularizer. The default value
%   is 1.
%   beta: rate parameter of exponential distribution. The default value
%   is 1.
%   tau: threshold parameter for obtaining functional units. The default
%   value is 0.3.
%   tol: tolerance of iteration. The default value is 0.01.
%   n_iter: the number of iterations limited in RSRGM. The default value is
%   150.
%   F_init: the initial value of protein-complex indication matrix F.

% Outputs:
%   F: protein-complex indication matrix.
%   F_star: resultant protein-complex indication matrix.
%   s: value of the objective function of (8).

    if nargin < 8
        F_init = rand(size(W,1),K);
    end

    if nargin < 7
        n_iter = 150;
    end
    
    if nargin < 6
        tol = 0.01;
    end
    
    if nargin < 5
        tau = 0.3;
    end
    
    if nargin < 4
        beta = 1;
    end
    
    if nargin < 3
         lambda = 1; 
    end
    
    if nargin<2
        error('Yor need input adjacent matrix W and maximum possible number K of protein complexes.');
    end
    
    n = size(W,1);
    D = diag(sum(W,2));
    
    F_old = F_init; %Initialize matrix F randomly.
    
    
    for i  = 1: n_iter
        % Update F according to Equation (14).
        
        F = 0.5*F_old + 0.5*F_old.*(  ((W./(F_old*F_old'+eps))*F_old + lambda*W*F_old) ./ (ones(n,n)*F_old + 0.5*beta + lambda*D*F_old +eps) );   
        if norm(F - F_old) < tol
            break
        end
        F_old = F;
    end
    
    % Calculate the  value s of the objective function (8).  
    s = - sum( sum(W.*log(F*F'+eps) ) )  + sum(sum(F*F')) + beta*sum(sum(F)) + lambda*( trace(F'*D*F) - trace(F'*W*F) );
    
    % Obtain the resultant protein-group indication matrix F_star
    % according to Equation (17).
    F_star = F;
    F_star(F_star >= tau) =1;
    F_star(F_star < tau) =0;
    
    % Delete the columns of F_star that contain at most two nonzero
    % elements.
    F_star (:,sum(F_star)<= 2) = [];
    
    
    
