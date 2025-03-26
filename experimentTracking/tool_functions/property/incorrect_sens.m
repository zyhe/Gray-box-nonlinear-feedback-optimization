%% Analyze the invalidating effect of incorrect sensitivity
clear;
clc;
close all

% Set system matrices 
A = [0.5, 1; 0, 0.5;];
B = eye(2);
C = [1, -2; 0, 0.5];
H = C * ((eye(2) - A) \ B);
H_hat = [2, 0; 0, -1];

uInit = ones(2, 1);
xInit = zeros(2, 1);

obj = @(u) norm(H * u)^2/2;

% Implement the closed-loop response
itrCnt = 20;
eta = 0.1; % step size
objVals = zeros(2, itrCnt);

% use the true sensitivity
u = uInit;
x = xInit;
for i = 1:itrCnt
    objVals(1, i) = obj(u);

    % system dynamics
    y = C * x;
    x = A * x + B * u;
    
    % update inputs
    grad = H.' * y;
    u = u - eta * grad;
end

% use the inexact sensitivity
u = uInit;
x = xInit;
for i = 1:itrCnt
    objVals(2, i) = obj(u);

    % system dynamics
    y = C * x;
    x = A * x + B * u;
    
    % update inputs
    grad_inexact = H_hat.' * y;
    u = u - eta * grad_inexact;
end

%% plot the results
figure(1);
plot(1: itrCnt, objVals(1, :), 'LineWidth', 2.5);
hold on;
plot(1: itrCnt, objVals(2, :), 'LineWidth', 2.5, 'Color', '#4C0099');
hold off;
xlabel('Number of Iterations', 'fontsize', 16);
ylabel('Objective', 'fontsize', 16);
grid on;
grid minor;

fprintf('Finish implementing the closed-loop response.\n');
