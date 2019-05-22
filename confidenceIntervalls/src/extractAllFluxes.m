function [flux, error] = extractAllFluxes(lms, mass)
    flux = zeros(size(lms));
    error = zeros(size(lms));
    for i = 1:size(flux,1)
        for j = 1:size(flux,2)
            currentModel = lms{i,j};
            flux(i,j) = currentModel.Coefficients.Estimate(2)/mass;
            error(i,j) = currentModel.Coefficients.SE(2)/mass;
        end
    end
end

