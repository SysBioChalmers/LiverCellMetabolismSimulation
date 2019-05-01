function svals = interpolateMetabolite(T, X, S, tvals)
    mu = log(X(2)/X(1))/(T(2)-T(1));
    [mk] = polyfit(X,S,1);
    f = mu*mk(1);
    svals = S(1) + f/mu * X(1)*(exp(mu*(tvals'-T(1)))-1);
end

