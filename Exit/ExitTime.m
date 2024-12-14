function T = ExitTime(epsilon, x0, num_paths)
    if nargin < 3
        num_paths = 1000; % Default number of simulation paths
    end
    
    % Time discretization parameters
    dt = 0.001;
    T_max = 100; % Maximum simulation time
    N = ceil(T_max/dt);
    
    % Initialize arrays for storing paths and exit times
    exit_times = zeros(num_paths, 1);
    
    % Main simulation loop
    for path = 1:num_paths
        X = x0(1);
        Y = x0(2);
        t = 0;
        
        for i = 1:N
            % Calculate gradient of V
            [gradVx, gradVy] = gradV(X, Y);
            
            % Generate Wiener increments
            dW1 = sqrt(dt) * randn;
            dW2 = sqrt(dt) * randn;
            
            % Euler-Maruyama step
            X_new = X - gradVx * dt + sqrt(2*epsilon) * dW1;
            Y_new = Y - gradVy * dt + sqrt(2*epsilon) * dW2;
            
            % Check if process has reached origin (within tolerance)
            if norm([X_new, Y_new]) < 0.1
                exit_times(path) = t;
                break;
            end
            
            X = X_new;
            Y = Y_new;
            t = t + dt;
        end
        
        if t >= T_max
            exit_times(path) = T_max;
        end
    end
    
    % Calculate mean exit time
    T = mean(exit_times);
end

function [gradVx, gradVy] = gradV(x, y)
    % Calculate gradient of potential V
    I = eye(2);
    
    % Numerical gradient
    h = 1e-6;
    gradVx = -(log(mvnpdf([x+h, y], [1, 0], I) + mvnpdf([x+h, y], [-1, 0], I))/2 - ...
               log(mvnpdf([x-h, y], [1, 0], I) + mvnpdf([x-h, y], [-1, 0], I))/2)/(2*h);
    gradVy = -(log(mvnpdf([x, y+h], [1, 0], I) + mvnpdf([x, y+h], [-1, 0], I))/2 - ...
               log(mvnpdf([x, y-h], [1, 0], I) + mvnpdf([x, y-h], [-1, 0], I))/2)/(2*h);
end