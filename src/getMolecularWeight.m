function sumOfMass = getMolecularWeight(model, metabolite)
    metNames = modifyMetNames(model);
    massMap = getMassMap();
    curMet = findIndex(metNames, metabolite);
    [elements, useMat, exitFlag]=parseFormulas(model.metFormulas(curMet), true);
    sumOfMass = getMass(massMap, elements, useMat(1,:));
end

function sumOfMass = getMass(massMap, elements, useMat)
    elements = elements.abbrevs;
    sumOfMass = 0;
    for i = 1:length(useMat)
        if useMat(i) > 0
            mass = massMap(elements{i});
            sumOfMass = sumOfMass + useMat(i) * mass; 
        end
    end
end

function massMap = getMassMap()
    masses = [12.011
            14.00674
            15.9994
            32.066
            30.973762
            1.00794
            6.941
            22.989768
            24.305
            35.4527
            39.0983
            40.078
            55.845
            58.9332
            63.546
            65.39
            74.92159
            78.96
            79.904
            88.90585
            126.90447
            0
            0
            ];
    abr = {
    'C'
    'N'
    'O'
    'S'
    'P'
    'H'
    'Li'
    'Na'
    'Mg'
    'Cl'
    'K'
    'Ca'
    'Fe'
    'Co'
    'Cu'
    'Zn'
    'As'
    'Se'
    'Br'
    'Y'
    'I'
    'R'
    'X'};
    massMap = containers.Map(abr, masses);
end
