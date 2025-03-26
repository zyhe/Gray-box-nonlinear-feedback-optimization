clear;
clc;
close all

% Define the grid for important parameters
etaZerothOrderGrid = 1e-4; %[7.5e-5, 1e-4, 2.5e-4];
% etaHyrbidGrid = 2.5e-4; %[2.5e-4, 5e-4];
deltaGrid = 5e-2; %[1e-4, 1e-3, 1e-2, 1e-1];

% Fixed parameters
itrCnt = 8e4; % number of iterations
period = 1e3; % period of the time-varying disturbance
numInstance = ceil(itrCnt/period); % number of instances
numTrial = 30; % number of trials
etaFirstOrder = 5e-4; % step size for accurate models
etaFirstOrderIA = 5e-4; % step size for inaccurate models
etaFirstOrderRLS = 5e-4; % step size for sensitivity learning
coeff_Sigma_p = 1e-4; % coefficient for process noise covariance, bigger -> faster
coeff_Sigma_m = 5e-4; % coefficient for measurement noise covariance, smaller -> faster
sigma_u = 1e-4; % standard deviation of excitation noise

etaHybrid = 5e-4; % step size for gray-box pipeline
const = 5; % constant used in decaying weights; larger -> less accurate
etaES = 2.5e-5; % step size for extremum seeking
ap = 1e-6; % amplitude for extremum seeking   

% Generate parameter files
paramDir = 'params';
if ~exist(paramDir, 'dir')
    mkdir(paramDir);
end

paramIdx = 0;
for etaZerothOrder = etaZerothOrderGrid
    for delta = deltaGrid
        params = struct('itrCnt', itrCnt, 'period', period, 'numInstance', numInstance, 'numTrial', numTrial, ...
            'etaFirstOrder', etaFirstOrder, 'etaFirstOrderIA', etaFirstOrderIA, 'etaFirstOrderRLS', etaFirstOrderRLS, ...
            'coeff_Sigma_p', coeff_Sigma_p, 'coeff_Sigma_m', coeff_Sigma_m, 'sigma_u', sigma_u, ...
            'etaZerothOrder', etaZerothOrder, 'etaHybrid', etaHybrid, 'delta', delta, 'const', const, ...
            'etaES', etaES, 'ap', ap);
        
        paramFileName = fullfile(paramDir, sprintf('params_%d.mat', paramIdx));
        save(paramFileName, '-struct', 'params');
        paramIdx = paramIdx + 1;
    end
end

fprintf('%d Parameter files generated successfully.\n', ...
    length(etaZerothOrderGrid) * length(deltaGrid));
