function result = humanLiverAA(inputLabels)
labels =     {'Alanine'
    'Arginine'
    'Asparagine'
    'Aspartate'
    'Cysteine'
    'Glutamate'
    'Glutamine'
    'Glycine'
    'Histidine'
    'Isoleucine'
    'Leucine'
    'Lysine'
    'Methionine'
    'Phenylalanine'
    'Proline'
    'Serine'
    'Threonine'
    'Tryptophan'
    'Tyrosine'
    'Valine'};

values = [0.069
0.0427
0.0934*0.4
0.0934*0.6     
0.0182
0.1034*0.66 
0.1034*0.34
0.074
0.0249
0.0498
0.1023
0.0654
0.0163
0.0463
0.0518
0.0681
0.0575
0.0098
0.0242
0.0829
];

countLookUp = containers.Map(labels, values);

result = zeros(length(inputLabels),1);
for i = 1:length(inputLabels)
    result(i) = countLookUp(inputLabels{i});
end

end

