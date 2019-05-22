function degradation = estimateGlutamineDegradation(concentrations, rate)

deltaT = 1/(60*60);
degradation = zeros(length(concentrations),1);
for i = 1:length(degradation)
    degradation(i) = concentrations(i) * (exp(-rate*deltaT) - 1)/deltaT;
end


end

