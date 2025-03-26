function closedLoopSim(paramIdx)
    % Implement the closed-loop simulation involving a nonlinear system and different controllers
    % Input: id = 0 (default) -> run the simulation with the default parameters
    % Output: None
    
    clc;
    close all

    % Set the random seed
    seed = 15;
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

    % generate the problem instance
    problem = generateProblem(n, p, qy, qd);

    if nargin < 1
        % If paramIdx is not provided, set algorithmic parameters directly
        itrCnt = 4e4; % number of iterations
        numTrial = 30; % number of trials

        etaFirstOrder = 5e-4; % step size for accurate models
        etaFirstOrderIA = 5e-4; % step size for inaccurate models
        etaFirstOrderRLS = 5e-4; % step size for sensitivity learning
        coeff_Sigma_p = 1e-2; % coefficient for process noise covariance
        coeff_Sigma_m = 5e-4; % coefficient for measurement noise covariance
        sigma_u = 5e-5; % standard deviation of excitation noise

        etaZerothOrder = 1e-4; % step size for zeroth-order optimization
        etaHybrid = 2.5e-4; % step size for hybrid optimization
        delta = 3e-3; % smoothing parameter
        const = 5; % constant used in decaying weights; larger -> less accurate

        etaES = 5e-4; % step size for extremum seeking
        ap = 5e-2; % amplitude for extremum seeking   
    else
        % load algorithmic parameters based on the index
        paramFileName = fullfile(currentDir, 'params', sprintf('params_%d.mat', paramIdx));
        params = load(paramFileName);
        fprintf('Load parameters from %s\n', paramFileName);

        itrCnt = params.itrCnt;
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

        fprintf('The values of etaZeroth and delta are %f and %f\n', etaZerothOrder, delta);
        fprintf('Parameters loaded successfully.\n');
    end

    % model-based feedback optimization
    firstParams = struct('itrCnt', itrCnt, 'etaFirstOrder', etaFirstOrder, 'etaFirstOrderIA', etaFirstOrderIA, ...
            'etaFirstOrderRLS', etaFirstOrderRLS, 'coeff_Sigma_p', coeff_Sigma_p, 'coeff_Sigma_m', coeff_Sigma_m, 'sigma_u', sigma_u);

    % zeroth-order and gray-box feedback optimization
    zerothParams = struct('itrCnt', itrCnt, 'numTrial', numTrial, 'etaZerothOrder', etaZerothOrder, 'delta', delta);
    grayParams = struct('itrCnt', itrCnt, 'numTrial', numTrial, 'etaHybrid', etaHybrid, 'delta', delta, 'const', const);

    % extremum seeking
    esParams = struct('itrCnt', itrCnt, 'numTrial', numTrial, 'etaES', etaES, 'ap', ap);

    %% Closed-loop response
    fig = figure();
    
    % model-based FO
    [gradFirstOrder, gradFirstOrderIA, gradFirstRLS] = firstResponse(problem, firstParams, fig);

    % zeroth-order FO
    try
        [upperBound, lowerBound, avgResult] = zerothResponse(problem, zerothParams, fig);
    catch ME
        upperBound = [];
        lowerBound = [];
        avgResult = [];
        disp(['Error occurred for zeroth-order FO: ', ME.message]);
    end
    
    % gray-box FO
    try
        [upperBoundHybrid, lowerBoundHybrid, avgResultHybrid, alphaArr] = grayBoxResponse(problem, grayParams, fig);
    catch ME
        upperBoundHybrid = [];
        lowerBoundHybrid = [];
        avgResultHybrid = [];
        alphaArr = [];
        disp(['Error occurred for gray-box FO: ', ME.message]);
    end

    % extremum seeking
    [upperBoundES, lowerBoundES, avgResultES] = stochESResponse(problem, esParams, fig);


    %% Save results
    timeStamp = datestr(datetime('now'));
    timeStamp = regexprep(timeStamp,' ','-');
    timeStamp = regexprep(timeStamp,':','-');
    suffix = sprintf('-%.1e-%.1e-%.1e.mat', delta, etaZerothOrder, etaHybrid);
    fileName = [sprintf('job_%d-', paramIdx), timeStamp, suffix];

    save(fullfile(currentDir, 'data', fileName), ...
            'gradFirstOrder', 'gradFirstOrderIA', 'gradFirstRLS', ...
            'upperBound', 'lowerBound', 'avgResult', ...
            'upperBoundHybrid', 'lowerBoundHybrid', 'avgResultHybrid', 'alphaArr',...
            'upperBoundES', 'lowerBoundES', 'avgResultES'...
            );
    fprintf('Data saved to %s\n', fullfile(currentDir, 'data', fileName));

    % save(fullfile(currentDir, 'data', fileName), ...
    %         'upperBound', 'lowerBound', 'avgResult');

    saveas(fig, fullfile(currentDir, 'figure', [fileName, '.png']));
    fprintf('Figure saved to %s\n', fullfile(currentDir, 'figure', [fileName, '.png']));
    
    fprintf('Finish the program successfully.\n');
end
