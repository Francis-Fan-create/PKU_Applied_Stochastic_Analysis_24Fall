% (a) Calculate the specific heat of a q=3 Potts model using the Metropolis algorithm
% Parameters
N = 100;          % Lattice size NxN
q = 3;            % Number of states
J = 1;            % Coupling constant
kB = 1;           % Boltzmann constant
h = 0;            % External field
T_range = linspace(0.5, 1.5, 11);  % Changed: wider range around T* ≈ 0.995

% Initialize lattice randomly
lattice = randi(q, N, N);

% Main simulation loop
num_temps = length(T_range);
energy = zeros(1, num_temps);
specific_heat = zeros(1, num_temps);

for t = 1:num_temps
    beta = 1/(kB * T_range(t));
    
    % Increased equilibration and sampling
    n_equil = 100*N^2;     % Changed: more equilibration steps
    n_samples = 20*N^2;   % Changed: more sampling steps
    energy_samples = zeros(1, n_samples);
    

    % Metropolis steps
    for step = 1:(n_equil + n_samples)
        % Full lattice sweep
        i = randi(N);
        j = randi(N);
        new_spin = randi(q);
        
        delta_E = calculateEnergyDiff(lattice, i, j, new_spin, J, h);
        
        if delta_E <= 0 || rand < exp(-beta * delta_E)
            lattice(i,j) = new_spin;
        end
        
        % Sample after equilibration
        if step > n_equil
            current_energy = calculateTotalEnergy(lattice, J, h);
            energy_samples(step-n_equil) = current_energy;
        end
    end
    
    % Calculate observables
    energy(t) = mean(energy_samples)/N^2;
    specific_heat(t) = (kB * beta^2 * var(energy_samples))/N^2;
    
    % Progress indicator
    fprintf('Temperature %d/%d completed\n', t, num_temps);
end

% Plot and save energy/specific heat
figure('Position', [100 100 800 600]);
subplot(2,1,1);
plot(T_range, energy, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Temperature (T)', 'FontSize', 12);
ylabel('Energy per site (u)', 'FontSize', 12);
title('Internal Energy vs Temperature', 'FontSize', 14);
grid on;

subplot(2,1,2);
plot(T_range, specific_heat, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Temperature (T)', 'FontSize', 12);
ylabel('Specific Heat (c)', 'FontSize', 12);
title('Specific Heat vs Temperature', 'FontSize', 14);
grid on;

savefig('energy_specific_heat.fig');

%(b) Calculate the correlation length of a q=3 Potts model using the Wolff algorithm
% Add parameters for magnetization study
T_values = [0.5, 1.0, 1.5];  % Selected temperatures
h_range = linspace(-5, 5, 11);  % Range of external field
magnetization = zeros(length(T_values), length(h_range));

% Loop over temperatures
for t_idx = 1:length(T_values)
    T = T_values(t_idx);
    beta = 1/(kB * T);
    
    % Loop over external field values
    for h_idx = 1:length(h_range)
        h = h_range(h_idx);
        lattice = randi(q, N, N);  % Fresh random start
        mag_samples = [];
        
        % Metropolis steps
        for step = 1:120*N^2
            i = randi(N);
            j = randi(N);
            new_spin = randi(q);
            
            delta_E = calculateEnergyDiff(lattice, i, j, new_spin, J, h);
            
            if delta_E <= 0 || rand < exp(-beta * delta_E)
                lattice(i,j) = new_spin;
            end
            
            if step > 100*N^2  % Sample after equilibration
                current_mag = sum(lattice,"all")/N^2;
                mag_samples(end+1) = current_mag;
            end
        end
        
        magnetization(t_idx, h_idx) = mean(mag_samples);
    end
end


% Plot and save magnetization
figure('Position', [100 100 600 500]);
hold on;
for t_idx = 1:length(T_values)
    plot(h_range, magnetization(t_idx, :), '-o', ...
         'LineWidth', 1.5, ...
         'MarkerSize', 8, ...
         'DisplayName', sprintf('T = %.1f', T_values(t_idx)));
end
xlabel('External Field (h)', 'FontSize', 12);
ylabel('Magnetization (m)', 'FontSize', 12);
title('Magnetization vs External Field', 'FontSize', 14);
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;
hold off;

savefig('magnetization.fig');


% (c) Calculate the correlation length of a q=3 Potts model
correlation_length = calculateCorrelations(N,q,J,kB,h,T_range);
% Plot and save correlation length
figure('Position', [100 100 600 500]);
plot(T_range, correlation_length, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Temperature (T)', 'FontSize', 12);
ylabel('Correlation Length (ξ)', 'FontSize', 12);
title('Correlation Length vs Temperature', 'FontSize', 14);
grid on;

savefig('correlation_length.fig');


% (d) Calculate the scaling exponents of a q=3 Potts model
calculateScalingExponents(T_range, specific_heat, correlation_length);