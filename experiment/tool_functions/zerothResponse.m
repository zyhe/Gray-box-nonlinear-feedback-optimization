function [upperBound, lowerBound, avgResult] = zerothResponse(problem, algParams, fig)
    % Implement the closed-loop response when the zeroth-order feedback optimization is applied
    %
    % Inputs:
    %   problem  - Struct of problem data (struct)
    %   algParams - Struct of algorithmic parameters (struct)
    %   fig  - figure handle for plotting the results (figure)
    %
    % Outputs:
    %   upperBound  - upper bound on the convergence measures (1-by-itrCnt vector)
    %   lowerBound  - lower bound on the convergence measures (1-by-itrCnt vector)
    %   avgResult   - average convergence measures (1-by-itrCnt vector)
    %

    % profile on
    % Add the tool_functions folder to the MATLAB path using a relative path
    % addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'tool_functions'));

    %% initialization
    % extract fields from the problem struct and define important functions
    [A, B, B2, E, C, D, dx, dy, M1, M2, m3, rho, uOpt, valOpt, sensMat, ...
        sensMat2, dxMat, sensMatIA, yss, Phi, PhiTilde, gradPhiTilde] = extractProblem(problem);

    % extract problem parameters
    n = size(A, 1); % state
    p = size(B, 2); % input
    qy = size(C, 1); % output
    qd = size(E, 2); % disturbance

    % extract algorithmic parameters
    itrCnt = algParams.itrCnt; % number of iterations
    numTrial = algParams.numTrial; % number of trials
    etaZerothOrder = algParams.etaZerothOrder; % step size
    delta = algParams.delta; % smoothing parameter

    % curFilePath = matlab.desktop.editor.getActiveFilename;
    % pathData = [erase(curFilePath,'zerothOrderFO.m'),'data\'];
    % load([pathData,'problemData-20-Feb-2024-17-29-54.mat']);

    %% implementation of zeroth-order feedback optimization
    fprintf('Start running zeroth-order FO\n');

    % store the results
    gradZerothOrder = zeros(numTrial, itrCnt);

    for num = 1:numTrial
        wZeroth = zeros(p, 1); % set the initial input
        % generate standard normal random vectors
        v = randn(p, 1); v = 1 / norm(v) * v;
        uZerothFormer = wZeroth + delta * v;
        x = zeros(n, 1); % set the initial state
        % measure the initial output
        yFormer = C * (A * x + B * uZerothFormer + ...
            rho * B2 * (sin(uZerothFormer) + uZerothFormer.^2) + E * dx) + D * dy;
        PhiFormer = Phi(uZerothFormer, yFormer);
        
        % realization of the system dynamics and feedback controller
        for itr = 1:itrCnt
            gradZerothOrder(num, itr) = gradPhiTilde(wZeroth);
            % generate the current standard normal random vector
            v = randn(p, 1); v = 1 / norm(v) * v;
            uZerothCur = wZeroth + delta * v;

            % implement the system dynamics
            [x, yCur] = dynamics(x, uZerothCur, A, B, rho, B2, C, D, E, dx, dy);

            % calculate the difference of function values and construct gradient estimates
            PhiCur = Phi(uZerothCur, yCur);
            objDifference = PhiCur - PhiFormer;     
            % obtain the new control input
            wZeroth = wZeroth - etaZerothOrder * p * objDifference / delta * v;
            PhiFormer = PhiCur;
        end

        info = sprintf('Finish the %d round of trial', num);
        disp(info);
    end

    % discard the rows containing the NaN elements
    gradZerothOrder(any(isnan(gradZerothOrder), 2), :) = [];

    upperBound = max(gradZerothOrder);
    lowerBound = min(gradZerothOrder);
    avgResult = mean(gradZerothOrder);

    %% plot the results
    figure(fig);
    idx = [1, 100: 100: itrCnt];
    semilogy(idx, abs(upperBound(idx)), 'w');
    hold on;
    semilogy(idx, abs(lowerBound(idx)), 'w');
    hold on;
    fill([idx, fliplr(idx)], [abs(upperBound(idx)), fliplr(abs(lowerBound(idx)))],...
            [229 204 255] / 255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
    mf = semilogy(idx, abs(avgResult(idx)), 'LineWidth', 2.5, 'Color', '#4C0099');
    % legend([mf], {'Model-free'}, 'Location', 'northeast', 'fontsize', 14);
    % hold on;
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('$\|\nabla \tilde{\Phi}(u)\|^2$', 'Interpreter', 'latex', 'fontsize', 16);
    grid on;
    grid minor;

    % profile viewer

    %% save the generated data
    % suffix = sprintf('-%.3f-%.1e.mat', delta, etaZerothOrder);
    % timeStamp = datestr(datetime('now'));
    % timeStamp = regexprep(timeStamp, ' ', '-');
    % timeStamp = regexprep(timeStamp, ':', '-');
    % fileName = ['zerothFOdata-', timeStamp, suffix];
    % save([pathData, fileName], 'upperBound', 'lowerBound', 'avgResult');

    fprintf('Finish running zeroth-order FO\n');
end


% function y = prox(x, mu)
%     % implementation of the proximal operator associated with the scaled l1-norm function
%     y = sign(x) .* max(abs(x) - mu, 0);
% end
