function halfLife = loadHalflife(outputMets)
    halfLife = zeros(length(outputMets));
    gluta = findIndex(outputMets, 'glutamine[s]');
    amonia = findIndex(outputMets, 'NH3[s]');
    pyroglutamate = findIndex(outputMets, '5-oxoproline[s]');
    
    glutHalf = 0.0023;
    
    halfLife(gluta, gluta) = -glutHalf;
    halfLife(amonia, gluta) = glutHalf;
    halfLife(pyroglutamate, gluta) = glutHalf;
    
end

