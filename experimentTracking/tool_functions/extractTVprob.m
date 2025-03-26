function [A, B, B2, E, C, D, rho, dx, dy, quadMat, linearMat, sensMat, sensMat2, sensMatIA, dxMat, uOpt, valOpt, lb, ub, numInstance] = extractTVprob(problem)
    % Extract fields from the problem struct and store them into corresponding variables
    %
    % Inputs:
    %   problem - Struct of problem data (struct)
    %
    % Outputs:
    %   A, B, B2, E, C, D, rho, dx, dy, quadMat, linearMat, sensMat, sensMat2, dxMat, uOpt, valOpt, lb, ub, numInstance - Corresponding variables

    A = problem.A;
    B = problem.B;
    B2 = problem.B2;
    E = problem.E;
    C = problem.C;
    D = problem.D;
    rho = problem.rho;
    dx = problem.dx;
    dy = problem.dy;
    quadMat = problem.quadMat;
    linearMat = problem.linearMat;
    sensMat = problem.sensMat;
    sensMat2 = problem.sensMat2;
    sensMatIA = problem.sensMatIA;
    dxMat = problem.dxMat;
    uOpt = problem.uOpt;
    valOpt = problem.valOpt;
    lb = problem.lb;
    ub = problem.ub;
    numInstance = problem.numInstance;
end
