function [min_F, min_F_star, min_s] = multi_RSGNM(network, K, lambda, beta, tau, tol, run_times, n_iter , F_init)
% multi_RSGNM repeats the entire calculation of RSGNM multiple times and
% chooses the result that gives the lowest value of objective function of
% (8). 
%
% Inputs:
%   network: a structure that contains information of PPI network. Filed of
%   adjacent matrix is the adjacent matrix PPI network and filed of protein_list is the
%   name list of correspding proteins.
%   K: maximum number of possible protein complexes.
%   lambda: penalization parameter for Laplacian regularizer. The default value
%   is 1.
%   beta: rate parameter of exponential distribution. The default value
%   is 1.
%   tau: threshold parameter for obtaining functional units. The default
%   value is 0.3.
%   tol: tolerance of iteration. The default value is 0.01;
%   run_times: the number of times that we repreat the entire calculation
%   of RSGNM. The default value is 50.
%   n_iter: the number of iterations limited in RSGNM. The default value is
%   150.
%   F_init: the initial value of protein-complex indication matrix F. The defult value is a n by K matrix by random.
% 
% Outputs:
%   min_F: the optimal value of F that corresponds the result giving
%   the lowest value of objective function of (8).
%   min_F_star: the optimal value of F_star that corresponds the result giving
%   the lowest value of objective function of (8).
%   min_s: the optimal value of s that corresponds the result giving
%   the lowest value of objective function of (8).

    if nargin < 9
        F_init = rand(size(network.adjacent_matrix,1),K);
    end
    
    if nargin < 8
        n_iter = 150;
    end
    
    if nargin < 7
        run_times = 1;
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
        error('Yor need input the PPI network and maximum possible number K of protein complexes.');
    end
    
    min_s = inf;
    for i = 1:run_times
        
        disp(['This is the ' num2str(i) '-th run']);
        
        %Main algorithm of RSGNM described in Algorithm 1 
        [F, F_star, s] = RSGNM(network.adjacent_matrix, K, lambda, beta, tau, tol, n_iter, F_init);
        % Choose the result that gives the lowest value of objective
        % function of  (8).
        if s < min_s
            min_F = F;
            min_F_star = F_star;
            min_s = s;
        end
  
    end


    

    % Write the results of  protein complexes to file
    % 'cohesive_protein_complexes.txt', where each row corresponds to an identified  complex. 
    fid = fopen('CEARD_protein_complexes.txt','w');
    for i = 1:size(min_F_star,2)
        member_indices = find(min_F_star(:,i));
        for j = 1:length(member_indices)
            fprintf(fid, '%s\t', cell2mat(network.protein_list(member_indices(j))) );
        end
        fprintf(fid, '\n');
    end
    fclose(fid);