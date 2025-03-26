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
itrCnt = 3e4;


%% Plot results
f1 = figure();
set(f1,'position',[100 50 600 450]);

% first-order, model-based methods
figure(f1);
mb = semilogy(1: itrCnt, gradFirstOrder(1: itrCnt), 'LineWidth', 2.5, 'Color', '#0066CC');
hold on;
ix = semilogy(1: itrCnt, gradFirstOrderIA(1: itrCnt), 'LineWidth', 2.5, 'Color', '#A2142F');
hold on;
idx = [1, 100: 100: itrCnt];
sl = semilogy(idx, gradFirstRLS(idx), 'LineWidth', 2.5, 'Color', '#FF9933');
hold on;

% model-free method
semilogy(idx, abs(upperBound(idx)), 'w');
hold on;
semilogy(idx, abs(lowerBound(idx)), 'w');
hold on;
fill([idx, fliplr(idx)], [abs(upperBound(idx)), fliplr(abs(lowerBound(idx)))],...
        [229 204 255] / 255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
mf = semilogy(idx, abs(avgResult(idx)), 'LineWidth', 2.5, 'Color', '#4C0099');
hold on;

% gray-box method
semilogy(idx, abs(upperBoundHybrid(idx)), 'w');
hold on;
semilogy(idx, abs(lowerBoundHybrid(idx)), 'w');
hold on;
fill([idx, fliplr(idx)], [abs(upperBoundHybrid(idx)), fliplr(abs(lowerBoundHybrid(idx)))], ...
    [153 255 204] / 255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
gb = semilogy(idx, abs(avgResultHybrid(idx)), 'LineWidth', 2.5, 'Color', '#33cc80');  % green
hold on;

% extremum seeking
semilogy(idx, lowerBoundES(idx), 'w');
hold on;
semilogy(idx, upperBoundES(idx), 'w');
hold on;
fill([idx, fliplr(idx)], [upperBoundES(idx), fliplr(lowerBoundES(idx))],...
     [224 236 255]./255, 'FaceAlpha', 0.5, 'LineStyle', 'none');
es = semilogy(idx, avgResultES(idx), 'LineWidth', 2.5, 'Color', '#99CCFF');
hold on;

xlabel('Number of Iterations', 'fontsize', 16);
ylabel('$\|\nabla \tilde{\Phi}(u)\|^2$', 'Interpreter', 'latex', 'fontsize', 16);
ylim([1e-3, 1e2]);
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
