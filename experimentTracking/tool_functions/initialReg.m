function regInit = initialReg(u, sensMat, sensMat2, dxInit, dyInit, D, M1Init, M2Init, rho, valOptInit)
    % Calculate the initial dynamic regret
    yss = @(u) sensMat * u + rho * sensMat2 * (sin(u) + u.^2) + dxInit + D * dyInit;
    Phi = @(u,y) (u.') * M1Init * u + M2Init * u + norm(y)^2;
    PhiTilde = @(u) Phi(u, yss(u));
    regInit = abs(PhiTilde(u) - valOptInit);
end
