function S = getSestimate(data,volume)
    %assumes sorted by time
    S = zeros(length(data),1);

    %ignore first volume widdrawal since t=0
    volume(1,2) = volume(2,2);
    diffvolume = diff(volume);
    diffvolume = diffvolume(:,2);
    diffvolume = [diffvolume;0]; %last sample does not matter
    
    wells = unique(data(:,2));
    
    for i = 1:length(wells)
        wellIndx = find(data(:,2) == wells(i));
        cumulativeRemoval = 0;        
        for j=1:length(wellIndx)
            volume(j,2)
            S(wellIndx(j)) =  data(wellIndx(j),3) * volume(j,2) + cumulativeRemoval;
            sampleRemoval = diffvolume(j) * data(wellIndx(j),3);
            cumulativeRemoval = cumulativeRemoval -sampleRemoval;
        end
    end


end
