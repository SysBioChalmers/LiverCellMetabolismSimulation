function drawSquare(x,y, height, width, strVal, colorVal)
    uX = x + width;
    uY = y + height;
    xdata = [x x uX uX];
    ydata = [y uY uY y];
    patch(xdata, ydata, colorVal, 'LineStyle', 'none');
    ypos = (y + uY)/2;
    xpos = x + width*0.05;
    text(xpos,ypos,strVal, 'color', 'w');
end

