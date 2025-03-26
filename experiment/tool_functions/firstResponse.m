function [gradFirstOrder, gradFirstOrderIA, gradFirstRLS] = firstResponse(problem, algParams, fig)
    % Implement the closed-loop response when the first-order feedback optimization is applied
    %
    % Inputs:
    %   problem  - Struct of problem data (struct)
    %   algParams - Struct of algorithmic parameters (struct)
    %   fig  - figure handle for plotting the results (figure)
    %
    % Outputs:
    %   gradFirstOrder  - gradient for first-order FO with accurate models (1-by-itrCnt vector)
    %   gradFirstOrderIA - gradient for first-order FO with inaccurate models (1-by-itrCnt vector)
    %   gradFirstRLS   - gradient for first-order FO with sensitivity learning (1-by-itrCnt vector)
    %

    % initialization
    [A, B, B2, E, C, D, dx, dy, M1, M2, m3, rho, uOpt, valOpt, sensMat, ...
        sensMat2, dxMat, sensMatIA, yss, Phi, PhiTilde, gradPhiTilde] = extractProblem(problem);

    % extract problem parameters
    n = size(A, 1); % state
    p = size(B, 2); % input
    qy = size(C, 1); % output
    qd = size(E, 2); % disturbance

    % extract algorithmic parameters
    itrCnt = algParams.itrCnt; % number of iterations
    etaFirstOrder = algParams.etaFirstOrder; % step size for accurate models
    etaFirstOrderIA = algParams.etaFirstOrderIA; % step size for inaccurate models
    etaFirstOrderRLS = algParams.etaFirstOrderRLS; % step size for sensitivity learning
    coeff_Sigma_p = algParams.coeff_Sigma_p; % coefficient for process noise covariance; bigger -> faster
    coeff_Sigma_m = algParams.coeff_Sigma_m; % coefficient for measurement noise covariance; smaller -> faster
    sigma_u = algParams.sigma_u; % standard deviation of excitation noise

    % implementation of first-order feedback controller with accurate models
    fprintf('Start running first-order FO with accurate models\n');
    uFirst = zeros(p, 1); % set the initial input
    x = zeros(n, 1); % set the initial state
    gradFirstOrder = zeros(1, itrCnt);

    for itr = 1:itrCnt
        gradFirstOrder(itr) = gradPhiTilde(uFirst);
        [x, y] = dynamics(x, uFirst, A, B, rho, B2, C, D, E, dx, dy);
        % construct the true gradient
        grad = -3 * m3 * norm(uFirst) * uFirst + 2 * M1 * uFirst + (M2.') ...
                + 2 * (sensMat + rho * sensMat2 * diag(cos(uFirst) + 2 * uFirst)).' * y;
        uFirst = uFirst - etaFirstOrder * grad;
    end

    fprintf('Finish running first-order FO with accurate models\n');

    % implementation of first-order feedback controller with inaccurate models
    fprintf('Start running first-order FO with inaccurate models\n');
    uFirstIA = zeros(p, 1); % set the initial input
    x = zeros(n, 1); % set the initial state
    gradFirstOrderIA = zeros(1, itrCnt);

    for itr = 1:itrCnt
        gradFirstOrderIA(itr) = gradPhiTilde(uFirstIA);
        [x, y] = dynamics(x, uFirstIA, A, B, rho, B2, C, D, E, dx, dy);
        grad = -3 * m3 * norm(uFirstIA) * uFirstIA + 2 * M1 * uFirstIA + (M2.') + 2 * sensMatIA.' * y;
        uFirstIA = uFirstIA - etaFirstOrderIA * grad;
    end

    fprintf('Finish running first-order FO with inaccurate models\n');

    % implementation of first-order feedback controller with sensitivity learning
    fprintf('Start running first-order FO with sensitivity learning\n');
    uFormer = zeros(p, 1);
    uCur = uFormer;
    x = zeros(n, 1);
    yFormer = C * (A * x + B * uFormer + rho * B2 * (sin(uFormer) + uFormer.^2) + E * dx) + D * dy;
    gradFirstRLS = zeros(1, itrCnt);
    sensEstErr = zeros(1, itrCnt);
    sensMatEst = sensMatIA;
    hHat = sensMatEst(:);
    Sigma = 1e-5 * eye(qy * p);

    for itr = 1:itrCnt
        gradFirstRLS(itr) = gradPhiTilde(uCur);
        sensEstErr(itr) = norm(sensMatEst - sensMat - rho * sensMat2 * diag(cos(uCur) + 2 * uCur));

        % system dynamics
        [x, yCur] = dynamics(x, uCur, A, B, rho, B2, C, D, E, dx, dy);

        % sensitivity learning based on recursive least squares
        delta_u = uCur - uFormer;
        U_delta = kron(delta_u.', eye(qy));
        delta_y = yCur - yFormer;
        Sigma_p = coeff_Sigma_p * eye(p * qy);
        Sigma_m = coeff_Sigma_m * eye(qy);

        K = Sigma * (U_delta.') / (Sigma_m + U_delta * Sigma * (U_delta.'));
        res = delta_y - U_delta * hHat;
        innov = K * res;
        hHat = hHat + innov;
        Sigma = (eye(p * qy) - K * U_delta) * Sigma + Sigma_p;
        sensMatEst = reshape(hHat, [qy, p]);

        uFormer = uCur;
        yFormer = yCur;

        % update the input
        grad = -3 * m3 * norm(uCur) * uCur + 2 * M1 * uCur + (M2.') + 2 * sensMatEst.' * yCur;
        noise_u_pre = randn(p, 1);
        uCur = uCur - etaFirstOrderRLS * grad + sigma_u * noise_u_pre;
    end

    fprintf('Finish running first-order FO with sensitivity learning\n');

    % plot the results
    figure(fig);
    semilogy(1: itrCnt, gradFirstOrder, 'LineWidth', 2.5, 'Color', '#0066CC');
    hold on;
    semilogy(1: itrCnt, gradFirstOrderIA, 'LineWidth', 2.5, 'Color', '#A2142F');
    hold on;
    idx = [1, 100: 100: itrCnt];
    semilogy(idx, gradFirstRLS(idx), 'LineWidth', 2.5, 'Color', '#FFB366');
    hold on;
    % ylim([1e-5, 1e2]);
    % legend({'Model-based, $H_k$',...
    %     'Model-based, $\hat{H}$','Sensitivity Learn.'},...
    %     'Interpreter','latex','Location','northeast','fontsize',14);
    % hold on;
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('$\|\nabla \tilde{\Phi}(u)\|^2$', 'Interpreter', 'latex', 'fontsize', 16);
    ylim([1e-4, 1e2]);
    grid on;
    grid minor;

    figure();
    semilogy(1: itrCnt, sensEstErr, 'LineWidth', 2.5, 'Color', '#00994C');
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('Sensitivity Estimation Error', 'fontsize', 16);
    grid on;
    grid minor;

    %% save the generated data
    % timeStamp = datestr(datetime('now'));
    % timeStamp = regexprep(timeStamp, ' ', '-');
    % timeStamp = regexprep(timeStamp, ':', '-');
    % fileName = ['firstAlldata-', timeStamp, '.mat'];
    % save(fileName, 'gradFirstOrder', 'gradFirstOrderIA', 'gradFirstRLS');
end
