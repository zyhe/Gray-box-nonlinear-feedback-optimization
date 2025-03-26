clear;
clc;
close all

%% Preparation
% Get the current script's directory
currentDir = fileparts(mfilename('fullpath'));
pathData = fullfile(currentDir, 'data');

% Sort files by date
S = dir(fullfile(pathData, '*.mat'));
[~, idx] = sort([S.datenum], 'descend');
mostRecentFile = S(idx(1)).name;

% Load the most recent .mat file
load(fullfile(pathData, mostRecentFile));

% itrCnt = length(gradFirstOrder);
itrCnt = 8e4;


%% Plot results
f1 = figure();
set(f1,'position',[100 50 600 450]);

% first-order, model-based methods
figure(f1);
mb = semilogy(1: itrCnt, avgRegFirstOrder(1: itrCnt), 'LineWidth', 2.5, 'Color', '#0066CC');
hold on;
ix = semilogy(1: itrCnt, avgRegFirstOrderIA(1: itrCnt), 'LineWidth', 2.5, 'Color', '#A2142F');
hold on;
idx = [1, 100: 100: itrCnt];
sl = semilogy(idx, avgRegFirstOrderRLS(idx), 'LineWidth', 2.5, 'Color', '#FFB366');
hold on;


% model-free method
semilogy(idx, upperBoundAvgReg(idx), 'w');
hold on;
semilogy(idx, lowerBoundAvgReg(idx), 'w');
hold on;
fill([idx, fliplr(idx)], [upperBoundAvgReg(idx), fliplr(lowerBoundAvgReg(idx))],...
        [229 204 255]./255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
mf = semilogy(idx, avgRegretZeroth(idx), 'LineWidth', 2.5, 'Color', '#4C0099');
hold on;

% gray-box method
semilogy(idx, upperBoundAvgRegHb(idx), 'w');
hold on;
semilogy(idx, lowerBoundAvgRegHb(idx), 'w');
hold on;
fill([idx, fliplr(idx)], [upperBoundAvgRegHb(idx), fliplr(lowerBoundAvgRegHb(idx))],...
     [153 255 204] ./ 255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
gb = semilogy(idx, avgRegretHb(idx), 'LineWidth', 2.5, 'Color', '#33CC80');
hold on;

% extremum seeking
semilogy(idx, upperBoundAvgRegES(idx), 'w');
hold on;
semilogy(idx, lowerBoundAvgRegES(idx), 'w');
hold on;
fill([idx, fliplr(idx)], [upperBoundAvgRegES(idx), fliplr(lowerBoundAvgRegES(idx))],...
     [224 236 255]./255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
es = semilogy(idx, avgRegretES(idx), 'LineWidth', 2.5, 'Color', '#99CCFF');
hold on;

xlabel('Number of Iterations', 'fontsize', 16);
ylabel('$\textup{Reg}_T^d/T$', 'Interpreter', 'latex', 'fontsize', 16);
ylim([1e-2, 1e3]);
grid on;
grid minor;

legend([mb, ix, sl, es, mf, gb], ...
    {'Model-based, $H_k$', 'Model-based, $\hat{H}$', ...
     'SL', 'ES', 'Model-free', 'Gray-box'}, ...
    'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 12, ...
    'NumColumns', 2); % Set number of columns
set(gca, 'fontsize', 14, 'FontName', 'Times New Roman');


%% Save the figure
set(f1, 'PaperUnits', 'inches');
set(f1, 'PaperPosition', [0 0 6 4.5]); % [left, bottom, width, height]
saveas(f1, fullfile(currentDir, 'figure', 'result.pdf'));

fprintf('Figure saved to %s\n', fullfile(currentDir, 'figure', 'result.pdf'));
fprintf("Finish the program\n");
