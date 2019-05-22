function involvedRxns = makeMetMetMap(smallModel, smallSolution, mets)
    tresh = 10^-6;
    filetype = 'PDF';

    labels = modifyMetNames(smallModel);
    metMatch = ismember(labels, mets);
    labels= labels(metMatch);
    modifiedS = smallModel.S(metMatch,:);
    modifiedS = modifiedS .* repmat(smallSolution', length(labels),1);
    modifiedS = modifiedS * 1000;
    
    involvedRxns = sum(abs(modifiedS)>tresh)>=2;
    
    adjMatrix = zeros(length(labels));
    for i = 1:size(modifiedS,2)
       producers = find(modifiedS(:,i)<0);
       consumers = find(modifiedS(:,i)>0);
       for j = 1:length(consumers)
           c = consumers(j);
            for k = 1:length(producers)
                p = producers(k);
                adjMatrix(p,c) = adjMatrix(p,c)+modifiedS(c,i);
            end
       end
    end
        
    graph_to_dot('displayNetwork/test.dot', adjMatrix, labels);
    %Render network with python
    delete('displayNetwork/test.dot.pdf')
    filepath = [cd '/displayNetwork/test.dot'];
    filepath = strrep(filepath, '\', '/');
    
    commandStr = ['python displayNetwork/makePdf.py ' filepath ' ' filetype];
     [status, commandOut] = system(commandStr);
     commandOut
     if status==0
         fprintf('pass\n');
     else
         fprintf('fail\n');
     end
end



function graph_to_dot(filename, adj, labels)
    fid = fopen(filename, 'w');

    fprintf(fid, 'digraph G {\n');
%    fprintf(fid, 'rankdir=LR;\n');

%    fprintf(fid, 'center = 1;\n');
    fprintf(fid, 'size=\"%d,%d\";\n', 10, 10);


    for i = 1:length(labels)
        additional = ', fillcolor = "#5D9BD3", fontcolor = "#FFFFFF", style="filled,setlinewidth(0)"';
        fprintf(fid, '%d [ label = "%s", fontsize=25, shape="box" %s];\n', i, labels{i}, additional);
    end

%Put all edges

    for node1 = 1:length(adj)
        arcs = find(adj(node1,:));
        
        for node2 = arcs
            flux = adj(node1,node2);
            lineWidth = calculateLineWidth(flux);
            lineWidth = 3;
            arrowSize = 1;
            fprintf(fid, '%d -> %d [label = "%2.0f", arrowsize=%2.2f, fontsize=20, color="#A0A0A0", penwidth=%f];\n', node1, node2, flux, arrowSize, lineWidth);  
        end
    end    
    fprintf(fid, '}\n'); 
    fclose(fid);
    
end


function lineWidth = calculateLineWidth(flux)
    lineWidth = 5*log10(flux);
    if lineWidth<0
        lineWidth = 0;
    end
    lineWidth = lineWidth + 0.5;
end
