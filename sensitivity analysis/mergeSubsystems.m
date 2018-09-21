function model = mergeSubsystems(model, A, newSub)
    affectedA = ismember(model.subSystems, A);
    model.subSystems(affectedA) = {newSub};
end

