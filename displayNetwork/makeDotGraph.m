function makeDotGraph(smallModel, smallSolution, mets) 
    fluxTresh = 1/10000;
    fluxTresh = 0;
    [cMatrix, labels, rxnStart] = generateBiPartiteFromMets(smallModel, smallSolution, mets, fluxTresh);
    graph_to_dot('displayNetwork/test.dot', cMatrix,  labels, rxnStart);
    %Render network with python
    delete('displayNetwork/test.dot.pdf')
    filetype = 'PDF';

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



function graph_to_dot(filename, adj, labels, rxnStart)
    adj = 1000 * adj;

    fid = fopen(filename, 'w');

    fprintf(fid, 'digraph G {\n');
    fprintf(fid, 'rankdir=LR;\n');

    fprintf(fid, 'center = 1;\n');
    fprintf(fid, 'size=\"%d,%d\";\n', 10, 10);


    for i = 1:rxnStart
        additional = ', fillcolor = "#5D9BD3", fontcolor = "#FFFFFF", style="filled,setlinewidth(0)"';
        fprintf(fid, '%d [ label = "%s", fontsize=25, shape="box" %s];\n', i, labels{i}, additional);
    end

    for i = (rxnStart+1):length(adj)
        in = find(adj(i,:));
        out = find(adj(:,i));
        if length(in) > 0 && length(out)
            additional = 'fillcolor = "#E0E0E0", shape="box", style="filled,setlinewidth(0)"';
        else
            additional = 'fillcolor = "#808080", fontcolor = "#FFFFFF", shape="box", style="filled,setlinewidth(0)"';
        end
        
        fprintf(fid, '%d [label = "%s", %s];\n', i, labels{i}, additional);
        
    end
    
lineWidth = 2;    


%Put all edges

    for node1 = 1:length(adj)
        arcs = find(adj(node1,:));
        
        for node2 = arcs
            flux = adj(node1,node2);
            lineWidth = calculateLineWidth(flux);
            arrowSize = 3;
            fprintf(fid, '%d -> %d [label = "%2.2f", arrowsize=%2.2f, fontsize=20, color="#A0A0A0", penwidth=%f];\n', node1, node2, flux, arrowSize, lineWidth);  
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
