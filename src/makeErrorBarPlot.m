function xvals = makeErrorBarPlot(data, errorLB, errorUB, labels, color)
hold all
totalWidth = 0.5;
errorWidth = 1.5*totalWidth/size(data,2);
% backgroundColor = [0.95 0.95 0.95
%                    0.9  0.9  0.9];
minVal = min(min([data errorLB errorUB]));
maxVal = max(max([data errorLB errorUB]));

for i = 1:size(data,1)    
    xvals = i + linspace(-totalWidth/2,totalWidth/2,size(data,2));
    
    for j = 1:size(data,2)
        x = xvals(j);
        y = data(i,j);
        eLB = y-errorLB(i,j);
        eUB = y-errorUB(i,j);
        bar(x, y, errorWidth, 'LineStyle','none', 'FaceColor', color(j,:))
        errorbar(x, y, eLB, eUB, 'k', 'LineStyle','none', 'HandleVisibility','off');
    end
end
xticks(1:size(data,1))
set(gca,'xticklabel',labels)
ylim([minVal maxVal])
xlim([0 length(labels)]+0.5)
xtickangle(90)
set(gca,'FontSize', 15, 'FontWeight', 'bold');
end

