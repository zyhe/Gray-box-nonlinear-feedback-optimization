function [upperBoundAvgRegES, lowerBoundAvgRegES, avgRegretES] = stochESTracking(problem, algParams, fig)
    % Implement the closed-loop response when the stochastic extremum seeking algorithm is applied
    %
    % Inputs:
    %   problem  - Struct of problem data (struct)
    %   algParams - Struct of algorithmic parameters (struct)
    %   fig  - figure handle for plotting the results (figure)
    %
    % Outputs:
    %   upperBoundAvgRegES  - upper bound on the average regret (1-by-itrCnt vector)
    %   lowerBoundAvgRegES  - lower bound on the average regret (1-by-itrCnt vector)
    %   avgRegretES         - average dynamic regret (1-by-itrCnt vector)
    %

    % profile on

    %% initialization
    [A, B, B2, E, C, D, rho, dx, dy, quadMat, linearMat, ...
        sensMat, sensMat2, sensMatIA, dxMat, uOpt, valOpt, lb, ub, numInstance] = extractTVprob(problem);
    [n, p] = size(B);
    qy = size(C, 1);
    qd = size(D, 2);

    % extract algorithmic parameters
    itrCnt = algParams.itrCnt; % number of iterations
    period = algParams.period; % period of the time-varying disturbance
    numTrial = algParams.numTrial; % number of trials
    etaES = algParams.etaES; % step size
    ap = algParams.ap; % amplitude

    % uInit = randn(p, 1);
    uInit = 0.5 * ones(p, 1);
    % uInit = zeros(p, 1);
    regInit = initialReg(uInit, sensMat, sensMat2, dx(:,1), dy(:,1), D, quadMat(:,:,1), linearMat(:,:,1), rho, valOpt(1));
    distInit = norm(uInit - uOpt(:, 1));

    %% implementation of extremum seeking
    fprintf('Start running stochastic extremum seeking\n');

    % store the results
    regretES = [regInit * ones(numTrial, 1), zeros(numTrial, itrCnt - 1)];
    trackErrPtES = [distInit * ones(numTrial, 1), zeros(numTrial, itrCnt - 1)];

    % closed-loop response
    for num = 1:numTrial
        wES = uInit;
        x = zeros(n, 1); % set the initial state

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

            % add perturbation to obtain input
            vSin = sin(randn(p, 1));
            uES = wES + ap * vSin;

            % implement the system dynamics
            [x, y] = dynamics(x, uES, A, B, rho, B2, C, D, E, dx, dy);

            % calculate the objective value
            PhiVal = Phi(uES, y);
            wES = wES - etaES * PhiVal * vSin;
            wES = projection(wES, lb, ub);

            % store convergence measures
            regretES(num, itr) = regretES(num, itr - 1) + abs(PhiTilde(wES) - valOpt(:, dIndex));
            trackErrPtES(num, itr) = norm(wES - uOpt(:, dIndex));
        end

        info = sprintf('Finish the %d round of trial', num);
        disp(info);
    end

    % discard the rows containing the NaN elements
    trackErrPtES(any(isnan(trackErrPtES), 2), :) = [];
    regretES(any(isnan(regretES), 2), :) = [];

    % results corresponding to average regret
    avgRegESData = regretES ./ (1:itrCnt);
    upperBoundAvgRegES = max(avgRegESData);
    lowerBoundAvgRegES = min(avgRegESData);
    avgRegretES = mean(avgRegESData);

    %% plot the results
    figure(fig);
    idx = [1, 10: 10: itrCnt];
    semilogy(idx, upperBoundAvgRegES(idx), 'w');
    hold on;
    semilogy(idx, lowerBoundAvgRegES(idx), 'w');
    hold on;
    fill([idx, fliplr(idx)], [upperBoundAvgRegES(idx), fliplr(lowerBoundAvgRegES(idx))],...
         [224 236 255]./255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
    semilogy(idx, avgRegretES(idx), 'LineWidth', 2.5, 'Color', '#99CCFF');
    hold on;
    % xlabel('Number of Iterations', 'fontsize', 16);
    % ylabel('$\|\nabla \tilde{\Phi}(u)\|^2$', 'Interpreter', 'latex', 'fontsize', 16);
    % grid on;
    % grid minor;

    %profile viewer
    % % save the generated data
    % timeStamp = datestr(datetime('now'));
    % timeStamp = regexprep(timeStamp,' ','-');
    % timeStamp = regexprep(timeStamp,':','-');
    % fileName = ['ESdata-',timeStamp,'.mat'];
    % save([pathData,fileName],'objValES','gradES');

    fprintf('Finish running stochastic extremum seeking\n');
end