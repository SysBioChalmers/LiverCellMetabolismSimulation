function graph_to_dot(filename, adj, labels, rxnStart, exchangeMets, biomassMets)
    adj = 1000 * adj;

    fid = fopen(filename, 'w');

    fprintf(fid, 'digraph G {\n');


    fprintf(fid, 'center = 1;\n');
    fprintf(fid, 'size=\"%d,%d\";\n', 10, 10);


    for i = 1:rxnStart
        additional = ', fillcolor = "#E0E0E0", style="filled,setlinewidth(0)"';

        if biomassMets(i)
            additional = ', fillcolor = "#A09F9F", style="filled,setlinewidth(0)"';
        end
        if exchangeMets(i)
            additional = ', fillcolor = "#5D9BD3", fontcolor = "#FFFFFF", style="filled,setlinewidth(0)"';
        end
        fprintf(fid, '%d [ label = "%s", fontsize=30, shape="box", height=0.3 %s];\n', i, labels{i}, additional);
    end

    allColapsedEdges = [];
    
    for i = (rxnStart+1):length(adj)
        attributes = ', fontsize=3, shape=box, fixedsize=true, width=0.2, height=0.2, fillcolor = "#FFFFFF", style="filled,setlinewidth(0)"';
        %attributes = ', shape=point';
        colaps = checkIfColapsable(adj, i);
        if isempty(colaps)
            fprintf(fid, '%d [ label = "%s" %s];\n', i, labels{i}, attributes);
        else
            allColapsedEdges = [allColapsedEdges;colaps];
        end
    end

    for node1 = 1:length(adj)
        arcs = find(adj(node1,:));
        colapsData = allColapsedEdges(node1 == allColapsedEdges(:,1),:);

        for node2 = arcs
            flux = adj(node1,node2);
    %       lineWidth = max(log(flux), 0.2); 
            lineWidth = 0.5;    

            if ismember(node2, colapsData)
                conector = colapsData(node2 == colapsData(:,2),3);
                newLabel = sprintf('%2.2f\n%s',flux, labels{node2});
                fprintf(fid, '%d -> %d [arrowsize=0.5, fontsize=3, color="#505050", penwidth=%f, label="%s"];\n', node1, conector, lineWidth, newLabel);                                  
            elseif ismember(node1, allColapsedEdges(:,2)) == false
                fprintf(fid, '%d -> %d [arrowsize=0.5, fontsize=3, color="#505050", penwidth=%f, label="%2.2f"];\n', node1, node2, lineWidth, flux);  
            end
        end
    end
    fprintf(fid, '}'); 
    fclose(fid);
end

function colaps = checkIfColapsable(adj, rxn)
    colaps = [];
    in = find(adj(:,rxn));
    out = find(adj(rxn,:));
    if  length(in) == 1 && length(out) == 1
       if adj(in, rxn) == adj(rxn, out)
           colaps = [in, rxn, out];
       end
    end
end

