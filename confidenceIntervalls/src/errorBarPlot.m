function errorBarPlot(i,y,err,color)
hold all
bar(i,y,'edgecolor', 'none', 'facecolor', color)
errorbar(i,y,err, 'k','HandleVisibility','off');
end

