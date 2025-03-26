function closedLoopTracking(paramIdx)
    % Implement the closed-loop simulation involving a nonlinear system and different controllers
    % Input: id = 0 (default) -> run the simulation with the default parameters
    % Output: None
    
    clc;
    close all

    % Set the random seed
    seed = 21;
    rng(seed);

    % Get the current script's directory
    currentDir = fileparts(mfilename('fullpath'));

    % Add the 'tool_functions' folder to the MATLAB path
    addpath(fullfile(currentDir, 'tool_functions'));

    % set problem parameters of dimensions
    n = 30; % state
    p = 15; % input
    qy = 10; % output
    qd = 10; % disturbance

    if nargin < 1
        % If paramIdx is not provided, set algorithmic parameters directly
        itrCnt = 8e4; % number of iterations
        period = 1e3; % period of the time-varying disturbance
        numInstance = ceil(itrCnt/period); % number of instances
        numTrial = 3; % number of trials

        etaFirstOrder = 5e-4; % step size for accurate models
        etaFirstOrderIA = 5e-4; % step size for inaccurate models
        etaFirstOrderRLS = 5e-4; % step size for sensitivity learning
        coeff_Sigma_p = 1e-4; % coefficient for process noise covariance, bigger -> faster
        coeff_Sigma_m = 5e-4; % coefficient for measurement noise covariance, smaller -> faster
        sigma_u = 1e-4; % standard deviation of excitation noise

        etaZerothOrder = 1e-4; % step size for zeroth-order optimization
        etaHybrid = 5e-4; % step size for gray-box pipeline
        delta = 5e-2; % smoothing parameter
        const = 5; % constant used in decaying weights; larger -> less accurate

        etaES = 2.5e-5; % step size for extremum seeking
        ap = 1e-6; % amplitude for extremum seeking
    else
        % load algorithmic parameters based on the index
        paramFileName = fullfile(currentDir, 'params', sprintf('params_%d.mat', paramIdx));
        params = load(paramFileName);
        fprintf('Load parameters from %s\n', paramFileName);

        itrCnt = params.itrCnt;
        period = params.period;
        numInstance = params.numInstance;
        numTrial = params.numTrial;
        etaFirstOrder = params.etaFirstOrder;
        etaFirstOrderIA = params.etaFirstOrderIA;
        etaFirstOrderRLS = params.etaFirstOrderRLS;
        coeff_Sigma_p = params.coeff_Sigma_p;
        coeff_Sigma_m = params.coeff_Sigma_m;
        sigma_u = params.sigma_u;
        etaZerothOrder = params.etaZerothOrder;
        etaHybrid = params.etaHybrid;
        delta = params.delta;
        const = params.const;
        etaES = params.etaES;
        ap = params.ap;
    end

    % generate the problem instance
    problem = generateTVprob(n, p, qy, qd, numInstance);

    % model-based feedback optimization
    firstParams = struct('itrCnt', itrCnt, 'period', period, 'etaFirstOrder', etaFirstOrder, 'etaFirstOrderIA', etaFirstOrderIA, ...
            'etaFirstOrderRLS', etaFirstOrderRLS, 'coeff_Sigma_p', coeff_Sigma_p, 'coeff_Sigma_m', coeff_Sigma_m, 'sigma_u', sigma_u);

    % zeroth-order and gray-box feedback optimization
    zerothParams = struct('itrCnt', itrCnt, 'period', period, 'numTrial', numTrial, 'etaZerothOrder', etaZerothOrder, 'delta', delta);
    grayParams = struct('itrCnt', itrCnt, 'period', period, 'numTrial', numTrial, 'etaHybrid', etaHybrid, 'delta', delta, 'const', const);

    % extremum seeking
    esParams = struct('itrCnt', itrCnt, 'period', period, 'numTrial', numTrial, 'etaES', etaES, 'ap', ap);

    %% Closed-loop response
    fig = figure();
    
    % model-based FO
    [avgRegFirstOrder, avgRegFirstOrderIA, avgRegFirstOrderRLS] = firstTracking(problem, firstParams, fig);

    % zeroth-order FO
    try
        [upperBoundAvgReg, lowerBoundAvgReg, avgRegretZeroth] = zerothTracking(problem, zerothParams, fig);
    catch ME
        upperBoundAvgReg = [];
        lowerBoundAvgReg = [];
        avgRegretZeroth = [];
        disp(['Error occurred for zeroth-order FO: ', ME.message]);
    end
    
    % gray-box FO
    try
        [upperBoundAvgRegHb, lowerBoundAvgRegHb, avgRegretHb, alphaArr] = grayBoxTracking(problem, grayParams, fig);
    catch ME
        upperBoundAvgRegHb = [];
        lowerBoundAvgRegHb = [];
        avgRegretHb = [];
        alphaArr = [];
        disp(['Error occurred for gray-box FO: ', ME.message]);
    end

    % extremum seeking
    try
        [upperBoundAvgRegES, lowerBoundAvgRegES, avgRegretES] = stochESTracking(problem, esParams, fig);
    catch ME
        upperBoundAvgRegES = [];
        lowerBoundAvgRegES = [];
        avgRegretES = [];
        disp(['Error occurred for ES: ', ME.message]);
    end

    %% Save results
    timeStamp = datestr(datetime('now'));
    timeStamp = regexprep(timeStamp,' ','-');
    timeStamp = regexprep(timeStamp,':','-');
    suffix = sprintf('-%.1e-%.1e-%.1e.mat', delta, etaZerothOrder, etaHybrid);
    fileName = [sprintf('job_%d-', paramIdx), timeStamp, suffix];

    save(fullfile(currentDir, 'data', fileName), ...
        'avgRegFirstOrder', 'avgRegFirstOrderIA', 'avgRegFirstOrderRLS', ...
        'upperBoundAvgReg', 'lowerBoundAvgReg', 'avgRegretZeroth', ...
        'upperBoundAvgRegHb', 'lowerBoundAvgRegHb', 'avgRegretHb', 'alphaArr', ...
        'upperBoundAvgRegES', 'lowerBoundAvgRegES', 'avgRegretES');
    fprintf('Data saved to %s\n', fullfile(currentDir, 'data', fileName));

    saveas(fig, fullfile(currentDir, 'figure', [fileName, '.png']));
    fprintf('Figure saved to %s\n', fullfile(currentDir, 'figure', [fileName, '.png']));
    
    fprintf('Finish the program successfully.\n');
end
