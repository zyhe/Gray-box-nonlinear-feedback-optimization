function [upperBoundHybrid, lowerBoundHybrid, avgResultHybrid, alphaArr] = grayBoxResponse(problem, algParams, fig)
    % Implement the closed-loop response when the gray-box feedback optimization is applied
    %
    % Inputs:
    %   problem  - Struct of problem data (struct)
    %   algParams - Struct of algorithmic parameters (struct)
    %   fig  - figure handle for plotting the results (figure)
    %
    % Outputs:
    %   upperBoundHybrid  - upper bound on the convergence measures (1-by-itrCnt vector)
    %   lowerBoundHybrid  - lower bound on the convergence measures (1-by-itrCnt vector)
    %   avgResultHybrid   - average convergence measures (1-by-itrCnt vector)
    %   alphaArr          - the sequence of alpha coefficients (1-by-itrCnt vector)
    %

    % profile on
    % Add the tool_functions folder to the MATLAB path using a relative path
    % addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'tool_functions'));

    % Extract fields from the problem struct and define important functions
    [A, B, B2, E, C, D, dx, dy, M1, M2, m3, rho, uOpt, valOpt, sensMat, ...
        sensMat2, dxMat, sensMatIA, yss, Phi, PhiTilde, gradPhiTilde] = extractProblem(problem);

    % Extract problem parameters
    n = size(A, 1); % state
    p = size(B, 2); % input
    qy = size(C, 1); % output
    qd = size(E, 2); % disturbance

    % Extract algorithmic parameters
    itrCnt = algParams.itrCnt; % number of iterations
    numTrial = algParams.numTrial; % number of trials
    etaHybrid = algParams.etaHybrid; % step size
    delta = algParams.delta; % smoothing parameter
    const = algParams.const; % constant used in decaying weights; larger -> less accurate

    % Rule of setting alpha_k
    alpha = @(itr) min(1, const ./ (itr + 1) .^ (2 / 3));

    %% Implementation of gray-box feedback optimization
    fprintf('Start running gray-box FO\n');

    % Store the results
    gradHybrid = zeros(numTrial, itrCnt);

    for num = 1:numTrial
        w = zeros(p, 1); % set the initial input
        % Generate standard normal random vectors
        v = randn(p, 1); v = 1 / norm(v) * v;
        uFormer = w + delta * v;
        x = zeros(n, 1); % set the initial state
        yFormer = C * (A * x + B * uFormer ...
                + rho * B2 * (sin(uFormer) + uFormer.^2) + E * dx) + D * dy; % measure the initial output
        PhiFormer = Phi(uFormer, yFormer);

        for itr = 1:itrCnt
            gradHybrid(num, itr) = gradPhiTilde(w);
            % Generate the current standard normal random vector
            v = randn(p, 1); v = 1 / norm(v) * v;
            % Obtain the perturbed input
            uCur = w + delta * v;

            % Implement the system dynamics
            [x, yCur] = dynamics(x, uCur, A, B, rho, B2, C, D, E, dx, dy);

            % Calculate the difference of function values and construct gradient estimates
            PhiCur = Phi(uCur, yCur);
            objDifference = PhiCur - PhiFormer;
            gradEst = p * objDifference / delta * v;
            % Use first-order information to calculate approximate gradients
            gradFirstIA = -3 * m3 * norm(uCur) * uCur + 2 * M1 * uCur + (M2.') ...
                            + 2 * sensMatIA.' * yCur;
            % Calculate the hybrid gradient
            alphaCur = alpha(itr);
            gradComb = alphaCur * gradFirstIA + (1 - alphaCur) * gradEst;
            % Obtain the new control input
            w = w - etaHybrid * gradComb;
            PhiFormer = PhiCur;
        end

        info = sprintf('Finish the %d round of trial', num);
        disp(info);
    end

    % Discard the rows containing the NaN elements
    gradHybrid(any(isnan(gradHybrid), 2), :) = [];
    upperBoundHybrid = max(gradHybrid);
    lowerBoundHybrid = min(gradHybrid);
    avgResultHybrid = mean(gradHybrid);

    %% Plot the results
    figure(fig);
    idx = [1, 100: 100: itrCnt];
    semilogy(idx, abs(upperBoundHybrid(idx)), 'w');
    hold on;
    semilogy(idx, abs(lowerBoundHybrid(idx)), 'w');
    hold on;
    fill([idx, fliplr(idx)], [abs(upperBoundHybrid(idx)), fliplr(abs(lowerBoundHybrid(idx)))], ...
        [153 255 204] / 255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
    semilogy(idx, abs(avgResultHybrid(idx)), 'LineWidth', 2.5, 'Color', '#33cc80', 'DisplayName', 'Gray-Box');  % green
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('$\|\nabla \tilde{\Phi}(u)\|^2$', 'Interpreter', 'latex', 'fontsize', 16);
    grid on;
    grid minor;
    % legend([gb], {'Gray-Box'},'Location','northeast','fontsize',14);

    %% Analyze the coefficients
    figure();
    alphaArr = alpha(1:itrCnt);
    plot(1:itrCnt, alphaArr);
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('$\alpha$', 'Interpreter', 'latex', 'fontsize', 16);
    grid on;
    grid minor;

    % profile viewer

    %% Save the generated data
    % suffix = sprintf('-%.3f-%.1e.mat', delta, etaHybrid);
    % timeStamp = datestr(datetime('now'));
    % timeStamp = regexprep(timeStamp, ' ', '-');
    % timeStamp = regexprep(timeStamp, ':', '-');
    % fileName = ['hybridFOdata-', timeStamp, suffix];
    % save([pathData, fileName], 'upperBoundHybrid', 'lowerBoundHybrid', 'avgResultHybrid');

    fprintf('Finish running gray-box FO\n');
end


% function y = prox(x, mu)
%     % implementation of the proximal operator associated with the scaled l1-norm function
%     y = sign(x) .* max(abs(x) - mu, 0);
% end
