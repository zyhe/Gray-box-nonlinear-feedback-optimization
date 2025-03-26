function [x_next, y] = dynamics(x, u, A, B, rho, B2, C, D, E, dx, dy)
    % Dynamics of the nonlinear system
    %
    % Inputs:
    %   x  - Current state (n-by-1 vector)
    %   u  - Current input (p-by-1 vector)
    %   A  - State matrix (n-by-n matrix)
    %   B  - Input matrix (n-by-p matrix)
    %   rho  - Coefficient for the nonlinear term of the input (scalar)
    %   B2  - Input matrix for the nonlinear term (n-by-p matrix)
    %   C  - Output matrix (qy-by-n matrix)
    %   D  - Feedthrough matrix (qy-by-qd matrix)
    %   E  - Disturbance matrix (n-by-qd matrix)
    %   dx  - Disturbance (qy-by-1 vector)
    %   dy  - Feedthrough disturbance (qy-by-1 vector)
    %
    % Outputs:
    %   x_next  - Next state (n-by-1 vector)
    %   y  - Next output (qy-by-1 vector)
    %
    % x_next = A * x + B * u + rho * B2 * sin(u) + E * dx; % previous version
    x_next = A * x + B * u + rho * B2 * (sin(u) + u.^2) + E * dx;
    y = C * x_next + D * dy;
end