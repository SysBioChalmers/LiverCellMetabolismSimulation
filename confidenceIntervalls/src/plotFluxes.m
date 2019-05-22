function plotFluxes(flux,error,mets,conditions,color)
 k = 1;
 tickValues = zeros(length(mets),1);
    for i = 1:size(flux,1)
        tickValues(i) = k + size(flux,2)/2 -0.5;
        for j = 1:size(flux,2)
            errorBarPlot(k,flux(i,j),error(i,j),color(j,:))
            k = k+1;
        end

    end
    xticks(tickValues);
    xticklabels(mets);
    legend(conditions)
end

