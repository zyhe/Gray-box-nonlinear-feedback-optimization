function y = projection(x, lb, ub)
    % Projection onto the constraint set
    y = min(max(x, lb), ub);
    % y = x;
    % signIDlower = (x-lb)<0;
    % y(signIDlower) = lb(signIDlower);
    % signIDupper = (x-ub)>0;
    % y(signIDupper) = ub(signIDupper);
end