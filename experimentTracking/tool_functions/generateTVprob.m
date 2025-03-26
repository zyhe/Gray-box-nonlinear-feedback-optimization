function problem = generateTVprob(n, p, qy, qd, numInstance)
	% Generate system matrices and objectives for the time-varying optimization problem
    %
    % Inputs:
    %   n  - Dimension of the system (scalar, positive integer)
    %   p  - Dimension of the input (scalar, positive integer)
    %   qy - Dimension of the output (scalar, positive integer)
    %   qd - Dimension of the disturbance (scalar, positive integer)
	%   numInstance - Number of problem instances (scalar, positive integer)
    %
    % Outputs:
    %   problem  - Structure containing generated system matrices and optimization results
	%	


	% generate the system matrices
    [A, B, B2, E, C, D] = generateSystemMatrices(n, p, qy, qd);

	rho = 0; % coefficient before B2

	% generate disturbances
	[dx, dy] = generateDisturbances(qd, numInstance);

	% define the steady-state map
	[sensMat, sensMat2, dxMat] = defineSteadyStateMap(A, B, B2, E, C);
	sensMatIA = setInexactSensitivities(sensMat);

	% generate coefficient matrices in objectives
	[quadMat, linearMat] = generateCoefficientMatrices(p, numInstance);

	% specify the upper and lower bounds of constraint sets
	lb = -2 * rand(p,1);
	ub = 2 * rand(p,1);

	% obtain the trajectory of optimal solutions
	[uOpt, valOpt] = obtainOptimalSolutionsTV(n, p, numInstance, sensMat, sensMat2, dxMat, dx, dy, D, quadMat, linearMat, lb, ub, rho);

	% plot the trajectory of optimal solutions
	figure(1);
	plot(1: numInstance, vecnorm(uOpt));
	xlabel('Index of the instance');
	ylabel('Norm of the optimal solution');
	grid on
	grid minor

	figure(2);
	plot(1: numInstance, valOpt);
	xlabel('Index of the instance');
	ylabel('Optimal value');
	grid on
	grid minor

	%% save the generated data
	% curFilePath = matlab.desktop.editor.getActiveFilename;
	% pathData = [erase(curFilePath,'problemSetup.m'),'data\'];
	% timeStamp = datestr(datetime('now'));
	% timeStamp = regexprep(timeStamp,' ','-');
	% timeStamp = regexprep(timeStamp,':','-');
	% fileName = ['problemData-',timeStamp,'.mat'];
	% save([pathData,fileName],'A','B','B2','E','C','D','rho','dx','dy','quadMat','linearMat',...
	% 		'sensMat','sensMat2','dxMat','uOpt','valOpt','lb','ub'); 

	% store the generated data in a struct
	problem = struct('A', A, 'B', B, 'B2', B2, 'E', E, 'C', C, 'D', D, 'rho', rho, 'dx', dx, 'dy', dy, ...
			 'quadMat', quadMat, 'linearMat', linearMat, ...
			 'sensMat', sensMat, 'sensMat2', sensMat2, 'sensMatIA', sensMatIA, 'dxMat', dxMat, ...
			 'uOpt', uOpt, 'valOpt', valOpt, 'lb', lb, 'ub', ub, 'numInstance', numInstance);
end


function [A, B, B2, E, C, D] = generateSystemMatrices(n, p, qy, qd)
    unscaledA = randn(n);
    unscaledA = unscaledA - tril(unscaledA, -1) + triu(unscaledA, 1).'; % make a symmetric matrix
    A = 0.05 / max(eig(unscaledA)) * unscaledA; % adjust the spectral radius of A
    B = 0.5 * randn(n, p); % corresponds to the input
    B2 = 0.5 * randn(n, p); % corresponds to the sine of the input
    E = 0.5 * randn(n, qd); % corresponds to the disturbance
    C = 0.5 * randn(qy, n); % corresponds to the output
    D = 0.5 * randn(qy, qd); % corresponds to the disturbance
end

function [dx, dy] = generateDisturbances(qd, numInstance)
	dx = 2e-2 * randn(qd, numInstance);
	dy = 2e-2 * randn(qd, numInstance);
	% dx = 2e-3*rand(qd,numInstance) - 1e-3;
	% dy = 2e-3*rand(qd,numInstance) - 1e-3;
	% generate the time-varying variables
	% yTraj = 2e-1*randn(qy,numInstance);
end

function [sensMat, sensMat2, dxMat] = defineSteadyStateMap(A, B, B2, E, C)
    xssMat = eye(size(A, 1)) - A;
    sensMat = C * (xssMat \ B);
    sensMat2 = C * (xssMat \ B2);
    dxMat = C * (xssMat \ E);
end

function sensMatIA = setInexactSensitivities(sensMat)
	% set inexact sensitivities
	sensErr = 1.0 * mean(abs(sensMat), 'all');
	sensMatIA = sensMat + 2 * sensErr * rand(size(sensMat)) - sensErr;
end

function [quadMat, linearMat] = generateCoefficientMatrices(p, numInstance)
    quadMat = zeros(p, p, numInstance);
    linearMat = zeros(1, p, numInstance);

	% Main body of the coefficient matrices of the objective
	P = 0.5 * randn(p); 
	M1 = P * P.'; 
	M2 = 2 * ones(1, p) + 0.5 * randn(1, p);

	for t = 1:numInstance
		% generate coefficient matrices of the objective
		P = 1e-1 * randn(p); 
		perturb1 = P * P.'; 
		perturb2 = 2e-2 * randn(1, p);
		M1Cur = M1 + perturb1;
		M2Cur = M2 + perturb2;
		quadMat(:, :, t) = M1Cur;
		linearMat(:, :, t) = M2Cur;
	end
end

function [uOpt, valOpt] = obtainOptimalSolutionsTV(n, p, numInstance, sensMat, sensMat2, dxMat, dx, dy, D, quadMat, linearMat, lb, ub, rho);
	uOpt = zeros(p, numInstance);
	valOpt = zeros(1, numInstance);

	for t = 1:numInstance
		% access coefficient matrices
		M1Cur = quadMat(:, :, t);
		M2Cur = linearMat(:, :, t);

		% obtain the steady-state map
		yss = @(u) sensMat * u + rho * sensMat2 * (sin(u) + u.^2) + dxMat * dx(:, t) + D * dy(:, t);
		
		% define the objective function  
		Phi = @(u, y) (u.') * M1Cur * u + M2Cur * u + norm(y)^2;  % norm(y-yTraj(:,t))^2
		PhiTilde = @(u) Phi(u, yss(u));
		
		% obtain the optimal solution
		options = optimoptions('fmincon','display','off');
		% [uOptM,valOptM] = fmincon(PhiTilde,randn(p,1),[],[],[],[],lb,ub);
		[uOpt(:, t), valOpt(t)] = fmincon(PhiTilde, zeros(p, 1), [], [], [], [], lb, ub, [], options);

		% use YALMIP to find optimal solutions
		% yalmip('clear')
		% u = sdpvar(p,1);
		% Constraints = [lb <= u <= ub];
		% Objective = (u.')*M1*u + M2*u + norm(sensMat*u + rho*sensMat2*sin(u) + dxMat*dx(:,t) + D*dy(:,t))^2;
		% options = sdpsettings('verbose',1);
	    %% Solve the problem
		% sol = optimize(Constraints,Objective,options);
		% % Analyze error flags
		% if sol.problem == 0
		%  % Extract and display value
		%  uOptY = value(u);
	    %  valOptY = value(Objective);
		% else
		%  display('Hmm, something went wrong!');
		%  sol.info
		%  yalmiperror(sol.problem)
		% end
    end
    fprintf('Finish solving all the problem istances\n',t);
end