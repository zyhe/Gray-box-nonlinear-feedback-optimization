function[upperBoundAvgReg, lowerBoundAvgReg, avgRegretZeroth] = zerothTracking(problem, algParams, fig)
    % Implement the closed-loop response when the zeroth-order feedback optimization is applied
    %
    % Inputs:
    %   problem  - Struct of problem data (struct)
    %   algParams - Struct of algorithmic parameters (struct)
    %   fig  - figure handle for plotting the results (figure)
    %
    % Outputs:
    %   upperBoundAvgReg  - upper bound on the average regret (1-by-itrCnt vector)
    %   lowerBoundAvgReg  - lower bound on the average regret (1-by-itrCnt vector)
    %   avgRegretZeroth   - average dynamic regret (1-by-itrCnt vector)

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
    etaZerothOrder = algParams.etaZerothOrder; % step size
    delta = algParams.delta; % smoothing parameter

    uInit = 0.5 * ones(p, 1);
    % uInit = zeros(p, 1);
    regInit = initialReg(uInit, sensMat, sensMat2, dx(:,1), dy(:,1), D, quadMat(:,:,1), linearMat(:,:,1), rho, valOpt(1));
    distInit = norm(uInit - uOpt(:, 1));

    %% implementation of zeroth-order feedback optimization
    fprintf('Start running zeroth-order FO\n');
    regretZerothOrder = [regInit * ones(numTrial, 1), zeros(numTrial, itrCnt - 1)];
    trackErrPtZerothOrder = [distInit * ones(numTrial, 1), zeros(numTrial, itrCnt - 1)];

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

            % obtain the new control input
            w = w - etaZerothOrder * gradEst;
            w = projection(w, lb, ub);

            % store convergence measures
            trackErrPtZerothOrder(num, itr) = norm(w - uOpt(:, dIndex)); %/norm(uOpt(:,dIndex));
            regretZerothOrder(num, itr) = regretZerothOrder(num, itr-1) + abs(PhiTilde(w) - valOpt(dIndex));

            % prepare for the next iteration
            PhiFormer = PhiCur;
        end

        info = sprintf('Finish the %d round of trial',num);
        disp(info);
    end

    % discard the rows containing the NaN elements
    trackErrPtZerothOrder(any(isnan(trackErrPtZerothOrder), 2), :) = [];
    regretZerothOrder(any(isnan(regretZerothOrder), 2), :) = [];

    % results corresponding to average regret
    avgRegZerothData = regretZerothOrder ./ (1:itrCnt);
    upperBoundAvgReg = max(avgRegZerothData);
    lowerBoundAvgReg = min(avgRegZerothData);
    avgRegretZeroth = mean(avgRegZerothData);

    %% plot the results
    figure(fig);
    idx = [1, 10: 10: itrCnt];
    semilogy(idx, upperBoundAvgReg(idx), 'w');
    hold on;
    semilogy(idx, lowerBoundAvgReg(idx), 'w');
    hold on;
    fill([idx, fliplr(idx)], [upperBoundAvgReg(idx), fliplr(lowerBoundAvgReg(idx))],...
            [229 204 255]./255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
    semilogy(idx, avgRegretZeroth(idx), 'LineWidth', 2.5, 'Color', '#4C0099');
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('$\textup{Reg}_T^d/T$', 'Interpreter', 'latex', 'fontsize', 16);
    grid on;
    grid minor;

    fprintf('Finish running zeroth-order FO\n');
end
