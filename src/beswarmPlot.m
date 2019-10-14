function cordinates = beswarmPlot(data, radius)
    cordinates = cell(length(data),1);
    color = [0.8500    0.3250    0.0980];
    
    for i = 1:length(data)
        cordinates{i} = plotColumn(i, data{i}, radius);
        scatter(cordinates{i}, data{i}, 'filled', 'MarkerEdgeColor', color, 'MarkerFaceColor', color, 'MarkerFaceAlpha',.5);
        
        %scatter(data{i}, cordinates{i}, 'MarkerEdgeColor',color);
    end
    
    ylim([0 length(data)+1])
end


function cordinates = plotColumn(y, data, radius)
    cordinates = zeros(length(data),1);
    [data, indx] = sort(data);
    height = round(max(data)/radius)+1;
    width = length(data)+6;
    grid = zeros(height, width);
    leftLast = false;
    middle = round(width/2);
    
    for i = 1:length(data)
        xCord = data(i);
        xGrid = fix(xCord/radius)+1;
        yGrid = middle;
        distance = 0;

        while grid(xGrid, yGrid) == 1
            leftLast = not(leftLast);
            if leftLast
               yGrid = middle - ceil(distance);
            else
               yGrid = middle + ceil(distance);
            end
            distance = distance + 0.50001;
          
        end
        
        grid(xGrid, yGrid) = 1;
        yCord = y+(yGrid-middle)*radius;
        
        cordinates(i) = yCord;
    end

    %unsort data
    cordinates(indx) = cordinates;
end