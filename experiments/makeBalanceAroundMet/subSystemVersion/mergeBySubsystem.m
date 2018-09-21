function [uniqueIn, uniqueOut, newInflux, newOutflux] = mergeBySubsystem(model, inRxns, outRxns, influx, outflux)
eqnsIn = model.subSystems(inRxns);
eqnsOut = model.subSystems(outRxns);
uniqueIn = unique(eqnsIn, 'stable');
uniqueOut = unique(eqnsOut, 'stable');

newInflux = zeros(length(uniqueIn),1);
newOutflux = zeros(length(uniqueOut),1);

for i = 1:length(uniqueIn)
    affectedReactions = ismember(eqnsIn, uniqueIn{i});
    newInflux(i)=sum(influx(affectedReactions));
end

for i = 1:length(uniqueOut)
    affectedReactions = ismember(eqnsOut, uniqueOut{i});
    newOutflux(i)=sum(outflux(affectedReactions));
end

for i = 1:length(uniqueIn)
    outRxn = findIndex(uniqueOut, uniqueIn{i});
    if not(isempty(outRxn))
        In = newInflux(i);
        Out = newOutflux(outRxn);
        if In>Out
            newInflux(i) = In - Out;
            newOutflux(outRxn) = 0;
        else
            newOutflux(outRxn) = Out - In;
            newInflux(i) = 0;          
        end
    end
end


%remove empty
uniqueIn(newInflux == 0) = [];
newInflux(newInflux == 0) = [];
uniqueOut(newOutflux == 0) = [];
newOutflux(newOutflux == 0) = [];



%Merge small values
Treshold = 0.02;
Treshold = Treshold * sum(newInflux);
sumOfSmallIn = sum(newInflux(newInflux<Treshold));
sumOfSmallOut = sum(newOutflux(newOutflux<Treshold));

if sumOfSmallIn>10^-6
    uniqueIn(newInflux < Treshold) = [];
    newInflux(newInflux < Treshold) = [];
    uniqueIn{end+1} = 'Other';
    newInflux = [newInflux; sumOfSmallIn];
end

if sumOfSmallOut>10^-6
    uniqueOut(newOutflux < Treshold) = [];
    newOutflux(newOutflux < Treshold) = [];
    uniqueOut{end+1} = 'Other';
    newOutflux = [newOutflux; sumOfSmallOut];
end





end

