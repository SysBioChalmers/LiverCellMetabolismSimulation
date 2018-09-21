function graph_to_dot(filename, adj, labels, rxnStart, sourceMets, sinkMets)
    adj = 1000 * adj;

    fid = fopen([filename '.dot'], 'w');

    fprintf(fid, 'digraph G {\n');


    fprintf(fid, 'center = 1;\n');
    fprintf(fid, 'size=\"%d,%d\";\n', 10, 10);


    for i = 1:rxnStart
        additional = ', fillcolor = "#E0E0E0", style="filled,setlinewidth(0)"';

        if sinkMets(i)
            additional = ', fillcolor = "#909090", fontcolor = "#FFFFFF", style="filled,setlinewidth(0)"';
        end
        if sourceMets(i)
            additional = ', fillcolor = "#5D9BD3", fontcolor = "#FFFFFF", style="filled,setlinewidth(0)"';
        end
        fprintf(fid, '%d [ label = "%s", fontsize=25, shape="box" %s];\n', i, labels{i}, additional);
    end

    allColapsedEdges = [];

    
    
    for i = (rxnStart+1):length(adj)
        letter = char(64+i-rxnStart);
        %fprintf(fid, '%d [label = "", xlabel = "%d", shape=circle, fixedsize=true, width=0.1, height=0.1, fillcolor = "#A0A0A0", style="filled,setlinewidth(0)"];\n', i, i-rxnStart);
        fprintf(fid, '%d [label = "%s", shape=circle, fontcolor = "#FFFFFF", fontsize=15, fixedsize=true, width=0.3, height=0.3, fillcolor = "#A0A0A0", style="filled,setlinewidth(0)"];\n', i, letter);
        
    end
    
%Put all sources in a row
sourceIds = find(sourceMets);       
fprintf(fid, '{rank = same; ');
    for i = 1:length(sourceIds)
        fprintf(fid, '%i;', sourceIds(i));
    end
fprintf(fid, '}');
    
    
lineWidth = 1.8;    


%Put all edges

    for node1 = 1:length(adj)
        arcs = find(adj(node1,:));
        
        for node2 = arcs
            flux = adj(node1,node2);
            fprintf(fid, '%d -> %d [arrowsize=1, fontsize=8, color="#A0A0A0", penwidth=%f];\n', node1, node2, lineWidth);  
        end
    end    
    fprintf(fid, '}\n'); 
    fclose(fid);
    
    fid = fopen([filename '_legend.dot'], 'w');
    %Add legend
    nrOfRxns = length(adj)-rxnStart;
    fprintf(fid, 'digraph G {\n');  
    
       fprintf(fid, 'key [label=<<table border="0" cellpadding="2" cellspacing="0" cellborder="0" width="200px">\n');   
            for i = 1:nrOfRxns
                letter = char(64+i);
                fprintf(fid, '<tr><td align="left">%s. %s</td></tr>\n', letter, labels{rxnStart+i}); 
            end
			  
       fprintf(fid, '</table>>, shape="box", style="setlinewidth(0)"];\n');                
    
    fprintf(fid, '}\n');    
    fclose(fid);
end

