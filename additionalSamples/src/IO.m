function cellData = IO(filename)
    fidInFile = fopen(filename,'r');
    nextLine = fgets(fidInFile);
    cellData = {};
    i = 1;
    while nextLine >= 0
      nextLine = strtrim(nextLine);
      tmp = strsplit(nextLine,'\t', 'CollapseDelimiters', false);
      for j = 1:length(tmp)
        cellData{i,j} = tmp{j}; 
      end
      
      nextLine = fgets(fidInFile);    
      i = i + 1;
    end
    fclose(fidInFile);
end

