function [A, B, B2, E, C, D, dx, dy, M1, M2, m3, rho, uOpt, valOpt, sensMat, sensMat2, dxMat, sensMatIA, yss, Phi, PhiTilde, gradPhiTilde] = extractProblem(problem)
    % Extract fields from the problem struct and store them into corresponding variables
    %
    % Inputs:
    %   problem - Struct of problem data (struct)
    %
    % Outputs:
    %   A, B, B2, E, C, D, dx, dy, M1, M2, m3, rho, uOpt, valOpt, sensMat, sensMat2, dxMat, sensMatIA - Corresponding variables

    A = problem.A;
    B = problem.B;
    B2 = problem.B2;
    E = problem.E;
    C = problem.C;
    D = problem.D;
    dx = problem.dx;
    dy = problem.dy;
    M1 = problem.M1;
    M2 = problem.M2;
    m3 = problem.m3;
    rho = problem.rho;
    uOpt = problem.uOpt;
    valOpt = problem.valOpt;
    sensMat = problem.sensMat;
    sensMat2 = problem.sensMat2;
    dxMat = problem.dxMat;
    sensMatIA = problem.sensMatIA;

    % obtain the steady-state map
    yss = @(u) sensMat * u + rho * sensMat2 * (sin(u) + u.^2) + dxMat * dx + D * dy;

    % define the objective function
    Phi = @(u, y) -m3 * norm(u)^3 + (u.') * M1 * u + M2 * u + norm(y)^2;

    PhiTilde = @(u) Phi(u, yss(u));
    % second moment of the gradient
    gradPhiTilde = @(u) norm(-3 * m3 * norm(u) * u + 2 * M1 * u + M2.'...
                    + 2 * (sensMat + rho * sensMat2 * diag(cos(u) + 2 * u)).' * yss(u))^2;
end