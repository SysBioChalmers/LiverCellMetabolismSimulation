function makeErrorBarPlot(data, error, labels, color)
hold all
totalWidth = 0.5;
errorWidth = 1.5*totalWidth/size(data,2);
% backgroundColor = [0.95 0.95 0.95
%                    0.9  0.9  0.9];
minVal = min(min(data)) - max(max(error));
maxVal = max(max(data)) + max(max(error));


for i = 1:size(data,1)
%     bgColor = mod(i,2) + 1;
%     x = i + 0.5*[-1 1 1 -1];
%     y = [minVal minVal maxVal maxVal];
%     fill(x,y,backgroundColor(bgColor,:), 'LineStyle', 'none', 'HandleVisibility','off') 
    
    xvals = i + linspace(-totalWidth/2,totalWidth/2,size(data,2));
    for j = 1:size(data,2)
        x = xvals(j);
        y = data(i,j);
        e = error(i,j);
        bar(x, y, errorWidth, 'LineStyle','none', 'FaceColor', color(j,:))
        if e>0
            errorbar(x, y, e, 'k', 'LineStyle','none', 'HandleVisibility','off');
        end
    end
end
xticks(1:size(data,1))
set(gca,'xticklabel',labels)
ylim([minVal maxVal])
xlim([0 length(labels)]+0.5)
xtickangle(90)
set(gca,'FontSize', 15, 'FontWeight', 'bold');
end

