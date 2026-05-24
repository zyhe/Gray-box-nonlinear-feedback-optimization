clear;
clc;
close all

% Define the grid for important parameters
etaZerothOrderGrid = 1e-4; %[7.5e-5, 1e-4, 2.5e-4];
etaHyrbidGrid = 2.5e-4; %[2.5e-4, 5e-4];
deltaGrid = 3e-3; %[3e-3, 4e-3, 5e-3];

% Fixed parameters
itrCnt = 4e4;
numTrial = 40;
etaFirstOrder = 5e-4;
etaFirstOrderIA = 5e-4;
etaFirstOrderRLS = 5e-4;
coeff_Sigma_p = 1e-2;
coeff_Sigma_m = 5e-4;
sigma_u = 5e-5;
const = 1;
etaES = 5e-4;
ap = 5e-2;

% Generate parameter files
paramDir = 'params';
if ~exist(paramDir, 'dir')
    mkdir(paramDir);
end

paramIdx = 0;
for etaZerothOrder = etaZerothOrderGrid
    for etaHybrid = etaHyrbidGrid
        for delta = deltaGrid
            params = struct('itrCnt', itrCnt, 'numTrial', numTrial, ...
                'etaFirstOrder', etaFirstOrder, 'etaFirstOrderIA', etaFirstOrderIA, ...
                'etaFirstOrderRLS', etaFirstOrderRLS, 'coeff_Sigma_p', coeff_Sigma_p, ...
                'coeff_Sigma_m', coeff_Sigma_m, 'sigma_u', sigma_u, ...
                'etaZerothOrder', etaZerothOrder, 'etaHybrid', etaHybrid, ...
                'delta', delta, 'const', const, 'etaES', etaES, 'ap', ap);
            
            paramFileName = fullfile(paramDir, sprintf('params_%d.mat', paramIdx));
            save(paramFileName, '-struct', 'params');
            paramIdx = paramIdx + 1;
        end
    end
end

fprintf('%d Parameter files generated successfully.\n', ...
    length(etaZerothOrderGrid) * length(deltaGrid) * length(etaHyrbidGrid));
