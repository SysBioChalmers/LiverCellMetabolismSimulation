function color = getColor(i, total)
    maxValue = 0.65;
    minValue = 0.25;
    ratio = i/total;
    scaledRatio = ratio * (maxValue-minValue);
    colorFactor = minValue + scaledRatio;
    color = [colorFactor 0.1+colorFactor 0.85]; 
end

