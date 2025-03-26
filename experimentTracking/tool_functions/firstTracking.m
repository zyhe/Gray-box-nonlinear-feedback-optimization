function [avgRegFirstOrder, avgRegFirstOrderIA, avgRegFirstOrderRLS] = firstTracking(problem, algParams, fig)
    % Implement the closed-loop response when the first-order feedback optimization is applied
    %
    % Inputs:
    %   problem  - Struct of problem data (struct)
    %   algParams - Struct of algorithmic parameters (struct)
    %   fig  - figure handle for plotting the results (figure)
    %
    % Outputs:
    %   avgRegFirstOrder  - Average regret of the first-order FO with accurate models (1-by-itrCnt vector)
    %   avgRegFirstOrderIA - Average regret of the first-order FO with inexact models (1-by-itrCnt vector)
    %   avgRegFirstOrderRLS - Average regret of the first-order FO with recursive estimation (1-by-itrCnt vector)

    % initialization
    [A, B, B2, E, C, D, rho, dx, dy, quadMat, linearMat, ...
        sensMat, sensMat2, sensMatIA, dxMat, uOpt, valOpt, lb, ub, numInstance] = extractTVprob(problem);
    [n, p] = size(B);
    qy = size(C, 1);
    qd = size(D, 2);

    % extract algorithmic parameters
    itrCnt = algParams.itrCnt;
    period = algParams.period;
    etaFirstOrder = algParams.etaFirstOrder;
    etaFirstOrderIA = algParams.etaFirstOrderIA;
    etaFirstOrderRLS = algParams.etaFirstOrderRLS;
    coeff_Sigma_p = algParams.coeff_Sigma_p;
    coeff_Sigma_m = algParams.coeff_Sigma_m;
    sigma_u = algParams.sigma_u;

    % initial input
    uInit = 0.5 * ones(p, 1);
    % uInit = zeros(p, 1);
    regInit = initialReg(uInit, sensMat, sensMat2, dx(:,1), dy(:,1), D, quadMat(:,:,1), linearMat(:,:,1), rho, valOpt(1));
    distInit = norm(uInit - uOpt(:, 1));


    %% implementation of first-order feedback controller with accurate models
    fprintf('Start running first-order FO with accurate models\n');
    regretFirstOrder = [regInit, zeros(1, itrCnt - 1)];
    trackErrPtFirstOrder = [distInit, zeros(1, itrCnt - 1)];
    uFirst = uInit;
    x = zeros(n, 1);
    
    % realization of the system dynamics and feedback controller
    for itr = 2:itrCnt
        % cycle through the matrix containing the data of the disturbance
        dIndex = ceil(itr / period);
        dxCur = dx(:, dIndex);
        dyCur = dy(:, dIndex);
        M1 = quadMat(:, :, dIndex);
        M2 = linearMat(:, :, dIndex);
        yss = @(u) sensMat * u + rho * sensMat2 * (sin(u) + u.^2) + dxMat * dxCur + D * dyCur;
        % define the objective function
        Phi = @(u,y) (u.') * M1 * u + M2 * u + norm(y)^2;  % norm(y-yTraj(:,t))^2
        PhiTilde = @(u) Phi(u, yss(u));
        % implement dynamics
        [x, y] = dynamics(x, uFirst, A, B, rho, B2, C, D, E, dxCur, dyCur);
        % construct the gradient
        grad = 2 * M1 * uFirst + (M2.') + 2 * (sensMat + rho * sensMat2 * diag(cos(uFirst) + 2 * uFirst)).' * y;

        % iterative update
        uFirst = uFirst - etaFirstOrder * grad;
        uFirst = projection(uFirst, lb, ub);

        % store convergence measures
        regretFirstOrder(itr) = regretFirstOrder(itr - 1) + abs(PhiTilde(uFirst) - valOpt(:, dIndex));
        trackErrPtFirstOrder(itr) = norm(uFirst - uOpt(:, dIndex));
    end
    fprintf('Finish running first-order FO with accurate models\n');


    %% implementation of first-order feedback controller with inaccurate models
    fprintf('Start running first-order FO with inaccurate models\n');
    regretFirstOrderIA = [regInit, zeros(1, itrCnt - 1)];
    trackErrPtFirstOrderIA = [distInit, zeros(1, itrCnt - 1)];
    uFirstIA = uInit;
    x = zeros(n, 1);
    
    % realization of the system dynamics and feedback controller
    for itr = 2:itrCnt
        % cycle through the matrix containing the data of the disturbance
        dIndex = ceil(itr/period);
        dxCur = dx(:, dIndex);
        dyCur = dy(:, dIndex);
        M1 = quadMat(:, :, dIndex);
        M2 = linearMat(:, :, dIndex);
        yss = @(u) sensMat * u + rho * sensMat2 * (sin(u) + u.^2) + dxMat * dxCur + D * dyCur;
        % define the objective function
        Phi = @(u,y) (u.') * M1 * u + M2 * u + norm(y)^2;  % norm(y-yTraj(:,t))^2
        PhiTilde = @(u) Phi(u, yss(u));
        % implement dynamics
        [x, y] = dynamics(x, uFirstIA, A, B, rho, B2, C, D, E, dxCur, dyCur);
        % construct the inexact gradient
        grad = 2 * M1 * uFirstIA + (M2.') + 2 * (sensMatIA).' * y;

        % iterative update
        uFirstIA = uFirstIA - etaFirstOrderIA * grad;
        uFirstIA = projection(uFirstIA, lb, ub);

        % store convergence measures
        regretFirstOrderIA(itr) = regretFirstOrderIA(itr - 1) + abs(PhiTilde(uFirstIA) - valOpt(:, dIndex));
        trackErrPtFirstOrderIA(itr) = norm(uFirstIA - uOpt(:, dIndex));
    end
    fprintf('Finish running first-order FO with inaccurate models\n');

    
    %% implementation of first-order feedback controller with sensitivity learning
    fprintf('Start running first-order FO with sensitivity learning\n');
    regretFirstOrderRLS = [regInit, zeros(1, itrCnt - 1)];
    trackErrPtFirstOrderRLS = [distInit, zeros(1, itrCnt - 1)];
    uFormer = uInit;
    uCur = uFormer;
    x = zeros(n, 1);
    yFormer = zeros(qy, 1);
    sensMatEst = sensMatIA;
    hHat = sensMatEst(:); % qy*p-by-1 vector
    sensEstErr = zeros(1, itrCnt); % error of sensitivity learning
    Sigma = 1e-5 * eye(qy * p);

    % realization of the system dynamics and feedback controller
    for itr = 2:itrCnt
        % cycle through the matrix containing the data of the disturbance
        dIndex = ceil(itr/period);
        dxCur = dx(:, dIndex);
        dyCur = dy(:, dIndex);
        M1 = quadMat(:, :, dIndex);
        M2 = linearMat(:, :, dIndex);
        yss = @(u) sensMat * u + rho * sensMat2 * (sin(u) + u.^2) + dxMat * dxCur + D * dyCur;
        % define the objective function
        Phi = @(u,y) (u.') * M1 * u + M2 * u + norm(y)^2;  % norm(y-yTraj(:,t))^2
        PhiTilde = @(u) Phi(u, yss(u));
        % implement dynamics
        [x, yCur] = dynamics(x, uCur, A, B, rho, B2, C, D, E, dxCur, dyCur);
        
        % sensitivity learning based on recursive least squares
        delta_u = uCur - uFormer;
        U_delta = kron(delta_u.', eye(qy));
        delta_y = yCur - yFormer;
        Sigma_p = coeff_Sigma_p * eye(p * qy);
        Sigma_m = coeff_Sigma_m * eye(qy);

        K = Sigma * (U_delta.') / (Sigma_m + U_delta * Sigma * (U_delta.'));  % Kalman gain
        res = delta_y - U_delta * hHat;
        innov = K * res;
        hHat = hHat + innov;
        Sigma = (eye(p * qy) - K * U_delta) * Sigma + Sigma_p;
        sensMatEst = reshape(hHat, [qy, p]);

        uFormer = uCur;
        yFormer = yCur;

        % update the input
        grad = 2 * M1 * uCur + (M2.') + 2 * sensMatEst.' * yCur;
        noise_u_pre = randn(p, 1);
        uCur = uCur - etaFirstOrderRLS * grad + sigma_u * noise_u_pre;
        uCur = projection(uCur, lb, ub);

        % store convergence results
        regretFirstOrderRLS(itr) = regretFirstOrderRLS(itr - 1) + abs(PhiTilde(uCur) - valOpt(:, dIndex));
        trackErrPtFirstOrderRLS(itr) = norm(uCur - uOpt(:, dIndex));
        sensEstErr(itr) = norm(sensMatEst - sensMat - rho * sensMat2 * diag(cos(uCur) + 2 * uCur));
    end
    fprintf('Finish running first-order FO with sensitivity learning\n');

    avgRegFirstOrder = regretFirstOrder./(1:itrCnt);
    avgRegFirstOrderIA = regretFirstOrderIA./(1:itrCnt);
    avgRegFirstOrderRLS = regretFirstOrderRLS./(1:itrCnt);

    %% plot the results
    figure(fig);
    semilogy(1: itrCnt, avgRegFirstOrder, 'LineWidth', 2.5, 'Color', '#0066CC');
    hold on;
    semilogy(1: itrCnt, avgRegFirstOrderIA, 'LineWidth', 2.5, 'Color', '#A2142F');
    hold on;
    idx = [1, 100: 100: itrCnt];
    semilogy(idx, avgRegFirstOrderRLS(idx), 'LineWidth', 2.5, 'Color', '#FFB366');
    hold on;
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('$\textup{Reg}_T^d/T$', 'Interpreter', 'latex', 'fontsize', 16);
    % ylim([1e-3, 1e0]);
    grid on;
    grid minor;

    figure(4);
    semilogy(1: itrCnt, sensEstErr, 'LineWidth', 2.5, 'Color', '#00994C');
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('Sensitivity Estimation Error', 'fontsize', 16);
    grid on;
    grid minor;
end
