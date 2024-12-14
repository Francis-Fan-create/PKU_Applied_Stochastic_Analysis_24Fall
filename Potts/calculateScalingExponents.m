function [gamma, delta] = calculateScalingExponents(T_range, specific_heat, correlation_length)
    % Find T* (critical temperature) from max of specific heat
    [~, max_idx] = max(specific_heat);
    T_star = T_range(max_idx);
    
    % Calculate epsilon = |1 - T/T*|
    epsilon = abs(1 - T_range/T_star);
    
    % Convert to log space, removing zeros/infinities
    log_eps = log(epsilon+ 1e-10);
    log_c = log(specific_heat+ 1e-10);
    log_xi = log(correlation_length+ 1e-10);
    
    % Linear fits
    % Specific heat scaling (gamma)
    p_c = polyfit(log_eps, log_c, 1);
    gamma = -p_c(1);
    
    % Correlation length scaling (delta)
    p_xi = polyfit(log_eps, log_xi, 1);
    delta = -p_xi(1);
    
    % Plot results
    figure('Position', [100 100 800 800]);
    
    % Specific heat scaling
    subplot(2,1,1);
    scatter(log_eps, log_c, 'b', 'DisplayName', 'Data');
    hold on;
    plot(log_eps, polyval(p_c, log_eps), 'r--', 'DisplayName', 'Fit');
    xlabel('log(ε)');
    ylabel('log(c)');
    title(sprintf('Specific Heat Scaling: γ = %.2f', gamma));
    legend('show');
    grid on;
    
    % Correlation length scaling
    subplot(2,1,2);
    scatter(log_eps, log_xi, 'b', 'DisplayName', 'Data');
    hold on;
    plot(log_eps, polyval(p_xi, log_eps), 'r--', 'DisplayName', 'Fit');
    xlabel('log(ε)');
    ylabel('log(ξ)');
    title(sprintf('Correlation Length Scaling: δ = %.2f', delta));
    legend('show');
    grid on;

    savefig('scaling_exponents.fig');
end