% Clear workspace and close figures
clear all; clc; close all;
% Parameters for spatial investigation
epsilon = 0.5;  % Fixed epsilon
r = linspace(0.2, 2, 4);  % Radial distances
theta = linspace(0, 2*pi, 4);  % Angular positions
[R, THETA] = meshgrid(r, theta);
X = R.*cos(THETA);  % x-coordinates
Y = R.*sin(THETA);  % y-coordinates
num_paths = 500;    % Paths per point

% Initialize exit time matrix
T_values = zeros(size(X));

% Calculate exit times for different initial positions
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        x0 = [X(i,j), Y(i,j)];
        T_values(i,j) = ExitTime(epsilon, x0, num_paths);
        fprintf('Computing x0 = (%.2f, %.2f): T = %.4f\n', x0(1), x0(2), T_values(i,j));
    end
end

% Create visualization
figure;
[x_plot, y_plot] = meshgrid(linspace(-2,2,50), linspace(-2,2,50));
T_interp = griddata(X, Y, T_values, x_plot, y_plot, 'cubic');

% Plot contour map
contourf(x_plot, y_plot, T_interp, 20, 'LineColor', 'none');
colorbar;
hold on;
plot(0,0,'w.','MarkerSize',20); % Mark origin
axis equal;
xlabel('x');
ylabel('y');
title('Expected Exit Time T(ε,x_0) for ε = 0.5');
colormap('jet');

% Add polar grid for reference
for ri = r
    plot(ri*cos(theta), ri*sin(theta), 'w:', 'LineWidth', 0.5);
end
for thi = theta
    plot([0,2*cos(thi)], [0,2*sin(thi)], 'w:', 'LineWidth', 0.5);
end