% Clear workspace and close figures
clear all; clc; close all;
% Test the function with ε = 0.5 and initial point (1,0)
epsilon = 0.5;
x0 = [1, 0];
T = ExitTime(epsilon, x0);
fprintf('Expected exit time: %.4f\n', T);