function result = gemAA(inputLabels)
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

values = [0.0937
0.0507
0.04
0.04
0.0142
0.0613
0.0505
0.1831
0.0198
0.0309
0.0664
0.0571
0.0156
0.029
0.0853
0.0491
0.0402
0.0072
0.019
0.0471];

countLookUp = containers.Map(labels, values);

result = zeros(length(inputLabels),1);
for i = 1:length(inputLabels)
    result(i) = countLookUp(inputLabels{i});
end

end

