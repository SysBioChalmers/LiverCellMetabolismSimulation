function [xvals, yvals] = estimateExponential(data, N)
xvals = linspace(min(data(:,1)), max(data(:,1)), N);
yvals = zeros(N,1);
    for i = 1:(length(data)-1)
        X1 = data(i,1);
        X2 = data(i+1,1);
        Y1 = data(i,2);
        Y2 = data(i+1,2);        
        dt = X2-X1;
        mu = log(Y2/Y1)/dt
        xinterval = and(X1<=xvals, xvals<=X2);
        curX = xvals(xinterval)-X1;
        curY = Y1.*exp(mu.*curX);
        yvals(xinterval) = curY;
    end

end

