function [upperBoundES, lowerBoundES, avgResultES] = stochESResponse(problem, algParams, fig)
    % Implement the closed-loop response when the stochastic extremum seeking algorithm is applied
    %
    % Inputs:
    %   problem  - Struct of problem data (struct)
    %   algParams - Struct of algorithmic parameters (struct)
    %   fig  - figure handle for plotting the results (figure)
    %
    % Outputs:
    %   upperBoundES  - upper bound on the convergence measures (1-by-itrCnt vector)
    %   lowerBoundES  - lower bound on the convergence measures (1-by-itrCnt vector)
    %   avgResultES   - average convergence measures (1-by-itrCnt vector)

    % profile on

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
    etaES = algParams.etaES; % step size
    ap = algParams.ap; % amplitude

    %% implementation of extremum seeking
    fprintf('Start running stochastic extremum seeking\n');

    % store the results
    % objValES = zeros(1, itrCnt);
    gradES = zeros(1, itrCnt);

    % closed-loop response
    for num = 1:numTrial
        wES = zeros(p, 1);
        x = zeros(n, 1); % set the initial state

        for itr = 1:itrCnt
            gradES(num, itr) = gradPhiTilde(wES);
            % objValES(itr) = PhiTilde(wES);

            % add perturbation to obtain input
            vSin = sin(randn(p, 1));
            uES = wES + ap * vSin;

            % implement the system dynamics
            [x, y] = dynamics(x, uES, A, B, rho, B2, C, D, E, dx, dy);

            % calculate the objective value
            PhiVal = Phi(uES, y);
            wES = wES - etaES * PhiVal * vSin;
        end

        info = sprintf('Finish the %d round of trial', num);
        disp(info);
    end

    % convergence results
    lowerBoundES = min(gradES);
    upperBoundES = max(gradES);
    avgResultES = mean(gradES);

    %% plot the results
    % semilogy(1:itrCnt,abs(objValES-valOpt),'LineWidth',2.5,'Color','#0066CC');
    % hold on;
    figure(fig);
    idx = [1, 10: 10: itrCnt];
    semilogy(idx, lowerBoundES(idx), 'w');
    hold on;
    semilogy(idx, upperBoundES(idx), 'w');
    hold on;
    fill([idx, fliplr(idx)], [upperBoundES(idx), fliplr(lowerBoundES(idx))],...
         [224 236 255]./255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
    semilogy(idx, avgResultES(idx), 'LineWidth', 2.5, 'Color', '#99CCFF');
    hold on;
    xlabel('Number of Iterations', 'fontsize', 16);
    ylabel('$\|\nabla \tilde{\Phi}(u)\|^2$', 'Interpreter', 'latex', 'fontsize', 16);
    grid on;
    grid minor;

    %profile viewer
    % % save the generated data
    % timeStamp = datestr(datetime('now'));
    % timeStamp = regexprep(timeStamp,' ','-');
    % timeStamp = regexprep(timeStamp,':','-');
    % fileName = ['ESdata-',timeStamp,'.mat'];
    % save([pathData,fileName],'objValES','gradES');

    fprintf('Finish running stochastic extremum seeking\n');
end