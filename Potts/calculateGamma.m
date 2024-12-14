function Gamma_k = calculateGamma(lattices, k_start)
    [N, ~, ~] = size(lattices);  % N×N lattice, T time steps
    Gamma_k = zeros(1,N-k_start+1);
    
    
    for k = 1:N-k_start+1
        correlation_sum = 0;
        % Calculate C(i,j) for all pairs at distance k
        for i1 = 1:N
            for j1 = 1:N
                % Time series for current position
                spin_i = lattices(i1,j1,:);
                
                % Horizontal correlations
                j2_right = mod(j1 + k - 1, N) + 1;
                spin_j_right = lattices(i1,j2_right,:);
                
                correlation_sum = correlation_sum + ...
                    mean(spin_i .* spin_j_right) - ...
                    mean(spin_i) * mean(spin_j_right);
                
                % Vertical correlations
                i2_down = mod(i1 + k - 1, N) + 1;
                spin_j_down = lattices(i2_down,j1,:);
                
                correlation_sum = correlation_sum + ...
                    mean(spin_i .* spin_j_down) - ...
                    mean(spin_i) * mean(spin_j_down);
                
            end
        end
        fprintf("k = %d, correlation_sum = %f\n", k, correlation_sum);
        % We do not normalize since its just a proportionality constant
        Gamma_k(1,k) = correlation_sum;
    end
end