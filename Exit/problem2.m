% Clear workspace and close figures
clear all; clc; close all;

% Parameters
x0 = [1, 0];
epsilon_values = [2.0, 1.5, 1.0, 0.75, 0.5, 0.25, 0.1];
num_paths = 1000;

% Initialize array for exit times
T_values = zeros(size(epsilon_values));

% Calculate exit times for different epsilon
for i = 1:length(epsilon_values)
    T_values(i) = ExitTime(epsilon_values(i), x0, num_paths);
    fprintf('ε = %.2f: Expected exit time = %.4f\n', epsilon_values(i), T_values(i));
end

% Plot results
figure;
semilogx(epsilon_values, T_values, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 15);
grid on;
xlabel('ε');
ylabel('Expected Exit Time T(ε,x_0)');
title('Expected Exit Time vs ε');

% Fit exponential curve to analyze relationship
fit_epsilon = linspace(min(epsilon_values), max(epsilon_values), 100);
p = polyfit(log(epsilon_values), log(T_values), 1);
fit_T = exp(p(2)) * fit_epsilon.^p(1);

hold on;
plot(fit_epsilon, fit_T, 'r--');
legend('Simulation Data', sprintf('Fitted curve: T ∝ ε^{%.2f}', p(1)));