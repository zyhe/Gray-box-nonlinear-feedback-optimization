clear;
clc;
close all

% investigate the trajectories of adaptive weights

T = 5e4; % number of iterations
x = 0:1:T;

adaWeight = @(x,a,b) a*(T-x).^b./(x.^b+a*(T-x).^b);
decayWeight = @(x,C) min(1, C./(x+1).^(2/3));

% the design used in the version of Feb. 2024
figure(3);
for C = 1e1:40:3e2
    plot(x, decayWeight(x,C));
    hold on
end

% as a increases, the curve moves from left to right
% as b increases, the curve becomes more steep
figure(1);
for a = 0.5:0.1:1
    for b = 2:2:6
        plot(x, adaWeight(x,a,b));
        hold on
    end
end

% the design used in the version of Oct. 2022
figure(2);
a = 1;
b = 2;
plot(x, adaWeight(x,a,b));
