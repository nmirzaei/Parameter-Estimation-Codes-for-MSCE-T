function g = prior_guesser(parNum,lb,ub,NumSample,coreID,seeder)
    guess = zeros(NumSample,parNum);
    % Generate a unique seed using coreID and current time
    if seeder==1
        seed = str2double(coreID) + sum(clock);
        rng(seed, 'twister');  % Set the random seed
    end
    % Define parameter bounds as a matrix [min, max] for each parameter
    param_bounds = [lb ub];      
    % Number of parameters
    n_params = size(param_bounds, 1);
    
    % Generate Latin Hypercube Samples in [0, 1]
    lhs_samples = lhsdesign(NumSample, parNum);
    
    % Scale the samples to the defined parameter bounds
    for i = 1:parNum
        min_val = param_bounds(i, 1);
        max_val = param_bounds(i, 2);
        guess(:, i) = lhs_samples(:, i) * (max_val - min_val) + min_val;
    end
    g = guess';
end
