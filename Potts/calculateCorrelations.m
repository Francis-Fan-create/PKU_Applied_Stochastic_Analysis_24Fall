function [correlation_length] = calculateCorrelations(N, q, J, kB, h, T_range)
    % Parameters
    k_start = floor(N/2)+1;
    num_samples = 20*N^2;
    num_equilibration = 100*N^2;
    correlation_length = zeros(size(T_range));
    
    % Loop over temperatures
    for t_idx = 1:length(T_range)
        T = T_range(t_idx);
        beta = 1/(kB * T);
        
        % Initialize lattice and storage
        lattice = randi(q, N, N);
        lattice_samples = zeros(N, N, num_samples);
        
        % Equilibration and sampling
        for step = num_equilibration + num_samples
            % Metropolis updates
            i = randi(N);
            j = randi(N);
            new_spin = randi(q);
            delta_E = calculateEnergyDiff(lattice, i, j, new_spin, J, h);
            if delta_E <= 0 || rand < exp(-beta * delta_E)
                lattice(i,j) = new_spin;
            end
            if step> num_equilibration
                lattice_samples(:,:,step-num_equilibration) = lattice;
            end
        end

        % Calculate correlation function
        Gamma_k = calculateGamma(lattice_samples, k_start);
        
        % Fit correlation length
        k_values = k_start:N;
        log_Gamma = log(abs(Gamma_k)+1e-10);
        p = polyfit(k_values, log_Gamma, 1);
        correlation_length(t_idx) = -1/p(1);
    end
end