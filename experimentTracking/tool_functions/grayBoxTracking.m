function[upperBoundAvgRegHb, lowerBoundAvgRegHb, avgRegretHb, alphaArr] = grayBoxTracking(problem, algParams, fig)
    % Implement the closed-loop response when the gray-box feedback optimization is applied
    %
    % Inputs:
    %   problem  - Struct of problem data (struct)
    %   algParams - Struct of algorithmic parameters (struct)
    %   fig  - figure handle for plotting the results (figure)
    %
    % Outputs:
    %   upperBoundAvgRegHb  - upper bound on the average regret (1-by-itrCnt vector)
    %   lowerBoundAvgRegHb  - lower bound on the average regret (1-by-itrCnt vector)
    %   avgRegretHb         - average dynamic regret (1-by-itrCnt vector)
    %   alphaArr            - the sequence of alpha coefficients (1-by-itrCnt vector)

    % initialization
    [A, B, B2, E, C, D, rho, dx, dy, quadMat, linearMat, ...
        sensMat, sensMat2, sensMatIA, dxMat, uOpt, valOpt, lb, ub, numInstance] = extractTVprob(problem);
    [n, p] = size(B);
    qy = size(C, 1);
    qd = size(D, 2);

    % extract algorithmic parameters
    itrCnt = algParams.itrCnt; % number of iterations
    period = algParams.period; % period of the time-varying disturbance
    numTrial = algParams.numTrial; % number of trials
    etaHybrid = algParams.etaHybrid; % step size
    const = algParams.const; % constant used in decaying weights; larger -> less accurate
    delta = algParams.delta; % smoothing parameter

    % rule of setting alpha_k
    alpha = @(itr) min(1, const ./ (itr + 1).^(1 / 6));
    % alpha = @(itr) (itrCnt-itr).^4./(itr.^4 + (itrCnt-itr).^4);
    % alpha = @(itr) 0*itr;

    uInit = 0.5 * ones(p, 1);
    % uInit = zeros(p, 1);
    regInit = initialReg(uInit, sensMat, sensMat2, dx(:,1), dy(:,1), D, quadMat(:,:,1), linearMat(:,:,1), rho, valOpt(1));
    distInit = norm(uInit - uOpt(:, 1));

    %% implementation of gray-box feedback optimization
    fprintf('Start running gray-box FO\n');
    regretHybrid = [regInit * ones(numTrial, 1), zeros(numTrial, itrCnt - 1)];
    trackErrPtHybrid = [distInit * ones(numTrial, 1), zeros(numTrial, itrCnt - 1)];

    for num = 1:numTrial
        w = uInit;
        % generate standard normal random vectors
        x = zeros(n, 1); % set the initial state
        PhiFormer = 0;
 
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
            
            % generate the current standard normal random vector
            v = randn(p, 1); v = 1 / norm(v) * v;
            uCur = w + delta * v;
            
            % implement dynamics
            [x, yCur] = dynamics(x, uCur, A, B, rho, B2, C, D, E, dxCur, dyCur);
            PhiCur = Phi(uCur, yCur);
            objDifference = PhiCur - PhiFormer; 
            gradEst = p * objDifference / delta * v;

            % use first-order information to calculate approximate gradients
            gradFirstIA = 2 * M1 * uCur + (M2.') + 2 * sensMatIA.' * yCur;
            % calculate the hybrid gradient
            alphaCur = alpha(itr);
            gradComb = alphaCur * gradFirstIA + (1 - alphaCur) * gradEst;

            % obtain the new control input
            w = w - etaHybrid * gradComb;
            w = projection(w, lb, ub);

            % store convergence measures
            trackErrPtHybrid(num, itr) = norm(w - uOpt(:, dIndex)); %/norm(uOpt(:,dIndex));
            regretHybrid(num, itr) = regretHybrid(num, itr-1) + abs(PhiTilde(w) - valOpt(dIndex));

            % prepare for the next iteration
            PhiFormer = PhiCur;
        end

        info = sprintf('Finish the %d round of trial',num);
        disp(info);
    end

    % discard the rows containing the NaN elements
    trackErrPtHybrid(any(isnan(trackErrPtHybrid), 2), :) = [];
    regretHybrid(any(isnan(regretHybrid), 2), :) = [];

    % results corresponding to average regret
    avgRegretHybrid = regretHybrid ./ (1:itrCnt);
    upperBoundAvgRegHb = max(avgRegretHybrid);
    lowerBoundAvgRegHb = min(avgRegretHybrid);
    avgRegretHb = mean(avgRegretHybrid);

    %% plot the results
    figure(fig);
    idx = [1, 10: 10: itrCnt];
    semilogy(idx, upperBoundAvgRegHb(idx), 'w');
    hold on;
    semilogy(idx, lowerBoundAvgRegHb(idx), 'w');
    hold on;
    fill([idx, fliplr(idx)], [upperBoundAvgRegHb(idx), fliplr(lowerBoundAvgRegHb(idx))],...
         [153 255 204] ./ 255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
    semilogy(idx, avgRegretHb(idx), 'LineWidth', 2.5, 'Color', '#33CC80');
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('$\textup{Reg}_T^d/T$', 'Interpreter', 'latex', 'fontsize', 16);
    grid on;
    grid minor;

    %% Analyze the coefficients
    figure();
    alphaArr = alpha(1: itrCnt);
    plot(1:itrCnt, alphaArr);
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('$\alpha$', 'Interpreter', 'latex', 'fontsize', 16);
    grid on;
    grid minor;

    fprintf('Finish running gray-box FO\n');
end
