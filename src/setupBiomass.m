function model = setupBiomass(model, GAM, maintainance)
    model = configureSMatrix(model, GAM, 'HumanGrowth', 'human_growthMaintainance[c]');
    model = configureSMatrix(model, 4.71, 'HumanGrowth', 'human_protein_pool[c]');
    model = configureSMatrix(model, 0.11, 'HumanGrowth', 'human_RNAPool[c]');
    model = configureSMatrix(model, 0.09, 'HumanGrowth', 'human_DNAPool[c]');
    model = configureSMatrix(model, 0.04, 'HumanGrowth', 'cholesterol[c]');
    model = configureSMatrix(model, 0.12, 'HumanGrowth', 'phosphatidylPool[c]');
    model = configureSMatrix(model, 0.05, 'HumanGrowth', 'fattyAcidPool[c]');
    model = configureSMatrix(model, 0.01, 'HumanGrowth', 'heparan sulfate[c]');
    model = configureSMatrix(model, 0.02, 'HumanGrowth', 'glycogen[c]');
    model = configureSMatrix(model, 0.557, 'HumanGrowth', 'metabolitePool[c]');
    
    metFilter = true;
    %Filter out low abundant metabolites
    if metFilter
        metTresh = 10^-4;
        metRxn = findIndex(model.rxns, 'metabolitePool');
        metFilter = abs(model.S(:,metRxn))<metTresh;
        model.S(metFilter,metRxn) = 0;
    end
    
    oddChainFA = false;
    if oddChainFA
        lactRxn = createRXNStuct(model, 'oddChainFA', '0.067821068 myristic acid[c] + 0.04473304 pentadecylic acid[c] + 0.715728716 palmitate[c] + 0.02020202 margaric acid[c] + 0.038961039 oleate[c] + 0.112554113 stearate[c] => oddFA[c]', -1000, 1000, 'Artificall reactions');
        model=addRxns(model,lactRxn,3,'c',true);
        model = configureSMatrix(model, 0, 'ApproxPhosphatidate', 'palmitate[c]');
        model = configureSMatrix(model, 2, 'ApproxPhosphatidate', 'oddFA[c]');
        
        model = configureSMatrix(model, 0, 'FattyAcidPool', 'margaric acid[c]');
        model = configureSMatrix(model, 0, 'FattyAcidPool', 'myristic acid[c]');
        model = configureSMatrix(model, 0, 'FattyAcidPool', 'oleate[c]');
        model = configureSMatrix(model, 0, 'FattyAcidPool', 'palmitate[c]');
        model = configureSMatrix(model, 0, 'FattyAcidPool', 'palmitolate[c]');
        model = configureSMatrix(model, 0, 'FattyAcidPool', 'pentadecylic acid[c]');
        model = configureSMatrix(model, 0, 'FattyAcidPool', 'stearate[c]');
        model = configureSMatrix(model, 3, 'FattyAcidPool', 'oddFA[c]');
        	
    end   
    
    model = setParam(model, 'lb', 'human_ATPMaintainance', maintainance);
end