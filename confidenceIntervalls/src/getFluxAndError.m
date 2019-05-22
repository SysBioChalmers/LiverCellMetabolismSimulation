function fluxAndError = getFluxAndError(lm, mass)
    flux = lm.Coefficients.Estimate(2);
    error = lm.Coefficients.SE(2);
    fluxAndError = [flux error]/mass;
end

