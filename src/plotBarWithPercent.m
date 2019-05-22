function plotBarWithPercent(results, eqnsOut, conditionName) 
exMap = [67 116 160
         80 137 188
         91 155 213
         151 185 224
         190 209 234]/255;
     
exMap = [67 116 160
         158 158 158
         91 155 213
         211 211 211
         190 209 234]/255;     
     
percentTresh = 100;     
     
colormap(exMap) 
bar(results',0.9,'stack', 'LineStyle', 'none');
l=legend(eqnsOut,'location', 'ne');
%l.PlotChildren = l.PlotChildren(length(l.PlotChildren):-1:1);

xticks(1:length(conditionName))
xticklabels(conditionName)

for i = 1:length(conditionName)
    percentages = results(:,i);
    percentages = 100*percentages/sum(percentages);
    xvalues = cumsum([0;results(1:(end-1),i)]);
    for j = 1:length(percentages)
        if percentages(j)>percentTresh
            text(i, xvalues(j)+0.1, sprintf('%2.0f%%', percentages(j)), 'FontSize', 12)
        end
    end
end
legend boxoff
set(gca,'FontSize',12) 
end

