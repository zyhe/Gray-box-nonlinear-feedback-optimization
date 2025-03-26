function problem = generateProblem(n, p, qy, qd)
    % Generates system matrices and objective functions.
    %
    % Inputs:
    %   n  - Dimension of the system (scalar, positive integer)
    %   p  - Dimension of the input (scalar, positive integer)
    %   qy - Dimension of the output (scalar, positive integer)
    %   qd - Dimension of the disturbance (scalar, positive integer)
    %
    % Outputs:
    %   p  - Structure containing generated system matrices and optimization results
    %        (A, B, B2, E, C, D, dx, dy, M1, M2, m3, rho, uOpt, valOpt, sensMat, sensMat2, dxMat, sensMatIA)
    %

    % generate the system matrices
    [A, B, B2, E, C, D] = generateSystemMatrices(n, p, qy, qd);

    % coefficient before B2
    rho = 0.1;

    % generate the disturbance
    dx = 0.1 * randn(qd, 1);
    dy = 0.1 * randn(qd, 1);

    % define the steady-state map
    [sensMat, sensMat2, dxMat, yss] = defineSteadyStateMap(A, B, B2, E, C, D, dx, dy, rho);

    % define the nonconvex objective function
    [Phi, PhiTilde, gradPhiTilde, M1, M2, m3] = defineObjectiveFunction(p, sensMat, sensMat2, dxMat, D, dy, rho, yss);

    % calculate optimal solutions
    [uOpt, valOpt] = fminunc(PhiTilde, zeros(p, 1));

    % check the quality of solutions
    initValGap = PhiTilde(zeros(p, 1)) - valOpt;
    initGradPhi = gradPhiTilde(zeros(p, 1));
    optValGap = PhiTilde(uOpt) - valOpt;
    optGradPhi = gradPhiTilde(uOpt);

    % set inexact sensitivities
    sensMatIA = setInexactSensitivities(sensMat);

    % % save the generated data
    % curFilePath = matlab.desktop.editor.getActiveFilename;
    % pathData = [erase(curFilePath,'problemSetup.m'),'data\'];
    % % save the generated data
    % timeStamp = datestr(datetime('now'));
    % timeStamp = regexprep(timeStamp,' ','-');
    % timeStamp = regexprep(timeStamp,':','-');
    % fileName = ['problemData-',timeStamp,'.mat'];
    % save([pathData,fileName],'A','B','B2','E','C','D','dx','dy','M1','M2','m3','rho',...
    %         'uOpt','valOpt','sensMat','sensMat2','dxMat','sensMatIA'); 

    % store the generated data in a struct
    problem = struct('A', A, 'B', B, 'B2', B2, 'E', E, 'C', C, 'D', D, ...
    'dx', dx, 'dy', dy, 'M1', M1, 'M2', M2, 'm3', m3, 'rho', rho, ...
    'uOpt', uOpt, 'valOpt', valOpt, 'sensMat', sensMat, 'sensMat2', sensMat2, ...
    'dxMat', dxMat, 'sensMatIA', sensMatIA);
end


function [A, B, B2, E, C, D] = generateSystemMatrices(n, p, qy, qd)
    unscaledA = randn(n);
    unscaledA = unscaledA - tril(unscaledA, -1) + triu(unscaledA, 1).'; % make a symmetric matrix
    A = 0.05 / max(eig(unscaledA)) * unscaledA; % adjust the spectral radius of A
    B = 0.5 * randn(n, p); % corresponds to the input
    B2 = 0.5 * randn(n, p); % corresponds to the sine of the input
    E = 0.5 * randn(n, qd); % corresponds to the disturbance
    C = 0.5 * randn(qy, n); % corresponds to the output
    D = 0.5 * randn(qy, qd); % corresponds to the disturbance
end

function [sensMat, sensMat2, dxMat, yss] = defineSteadyStateMap(A, B, B2, E, C, D, dx, dy, rho)
    xssMat = eye(size(A, 1)) - A;
    sensMat = C * (xssMat \ B);
    sensMat2 = C * (xssMat \ B2);
    dxMat = C * (xssMat \ E);
    yss = @(u) sensMat * u + rho * sensMat2 * (sin(u) + u.^2) + dxMat * dx + D * dy; % modify if we add a quadratic term
end

function [Phi, PhiTilde, gradPhiTilde, M1, M2, m3] = defineObjectiveFunction(p, sensMat, sensMat2, dxMat, D, dy, rho, yss)
    P = 0.5 * randn(p);
    M1 = P * P.';
    M2 = 0.5 * randn(1, p);
    m3 = 5e-2;
    Phi = @(u, y) -m3 * norm(u) ^ 3 + (u.') * M1 * u + M2 * u + norm(y) ^ 2;
    PhiTilde = @(u) Phi(u, yss(u));
    % second moment of gradients
    gradPhiTilde = @(u) norm(-3 * m3 * norm(u) * u + 2 * M1 * u + M2.' ...
        + 2 * (sensMat + rho * sensMat2 * diag(cos(u) + 2 * u)).' * yss(u))^2;
end

function sensMatIA = setInexactSensitivities(sensMat)
    sensErr = 0.15 * mean(sensMat, 'all');
    sensMatIA = sensMat + 2 * sensErr * rand(size(sensMat)) - sensErr;
end
