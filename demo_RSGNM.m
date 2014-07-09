% Choose a network from the following three networks to test RSGNM.

% database='Gavin';
% database='Krogan'
 database='Collins';


switch database
    
    case 'Gavin'
        load ./data/Gavin_network.mat
        K = 500;
        lambda = 4;
        beta = 8;
        tau = 0.3;
        tol = 0.01;
        run_times = 50;
        n_iter = 150;
        F_init = rand(size(Gavin_network.adjacent_matrix,1),K);
        [F, F_star, s] = multi_RSGNM(Gavin_network, K, lambda, beta, tau, tol, run_times, n_iter , F_init);
        
    case 'Krogan'
        load ./data/Krogan_network.mat
        K = 500;
        lambda = 8;
        beta = 8;
        tau = 0.3;
        tol = 0.01;
        run_times = 50;
        n_iter = 150;
        F_init = rand(size(Krogan_network.adjacent_matrix,1),K);
        [F, F_star, s] = multi_RSGNM(Krogan_network, K, lambda, beta, tau, tol, run_times, n_iter , F_init);
     
        
    case 'Collins'
        load ./data/Collins_network.mat
        K = 500;
        lambda = 4;
        beta = 8;
        tau = 0.3;
        tol = 0.01;
        run_times = 1;
        n_iter = 150;
        F_init = rand(size(Collins_network.adjacent_matrix,1),K);
        [F, F_star, s] = multi_RSGNM(Collins_network, K, lambda, beta, tau, tol, run_times, n_iter , F_init);
          
end


     