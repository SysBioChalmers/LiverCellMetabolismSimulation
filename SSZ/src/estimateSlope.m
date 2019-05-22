function [slope, intercept, SE] = estimateSlope(time, cellCount)
    mdlL = fitlm(time, cellCount);
    intercept = mdlL.Coefficients.Estimate(1);
    slope = mdlL.Coefficients.Estimate(2);
    SE = mdlL.Coefficients.SE;
end

