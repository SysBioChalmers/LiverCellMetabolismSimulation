function mLL = minFunction(beta, dataX, dataY, modelfun)
    params = length(beta)/3;
    I1 = 1:params;
    I2 = (params+1):(params*2);
    I3 = (2*params+1):(params*3);
    x0 = [beta(I1) beta(I2)];
    sig0 = beta(I3);
    
    predictedY = modelfun(x0, dataX);
    mLL = 0;
    for i = 1:length(predictedY)
        curData = dataX(i,2);
        sig = sig0(curData);
        e = dataY(i) - predictedY(i);
        mLL = mLL -(0.5 * log(2*pi)) -(0.5 * log(sig.^2)) -(0.5 * (e/sig).^2);
    end

    mLL = -mLL; %minimize instead of maximize
end