function [sequenceProt, sequences, labels] = loadSequenceData(file)
    sequenceFile = importdata(file);
    sequences = sequenceFile.data;
    sequenceProt = sequenceFile.textdata;
%labels ={'A' 'R' 'N' 'D' 'C' 'E' 'Q' 'G' 'H' 'I' 'L' 'K' 'M' 'F' 'P' 'S' 'T' 'W' 'Y' 'V'};


labels ={'Alanine'
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
    
end

