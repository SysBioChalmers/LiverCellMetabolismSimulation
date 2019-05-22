function result = getMolWeight(inputLabels)
values = [89.0935
    174.2017
    132.12
    133.1032
    121.159
    147.1299
    146.14
    75.0669
    155.1552
    131.1736
    131.1736
    146.1882
    149.2124
    165.19
    115.131
    105.093
    119.1197
    204.2262
    181.1894
    117.1469];

    labels = {
    'Alanine'
    'Arginine'
    'Asparagine'
    'Aspartate'
    'Cysteine'
    'Glutamine'
    'Glutamate'
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


    countLookUp = containers.Map(labels, values);

    result = zeros(length(inputLabels),1);
    for i = 1:length(inputLabels)
        result(i) = countLookUp(inputLabels{i});
    end

end