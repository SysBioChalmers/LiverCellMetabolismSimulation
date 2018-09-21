function linkMets(model, smallSolution)

for i = 1:length(model.rxns)
    subRxns = split(model.rxns{i},'+');
    model.rxns{i} = subRxns{1};
end

sources = {'glucose[c]', 'pyruvate', 'alanine', 'glutamine'};
sinks = {'alanine', 'L-lactate', 'aspartate', 'AKG', 'glutamate'};
banMets = {'AKG', 'Pi' '(R)-methylmalonyl-CoA', 'ubiquinol', 'GSH', 'ATP', 'AMP', 'ADP', 'CoA' 'H+', 'H2O', 'NADH', 'NAD+', 'GMP', 'CTP', 'PPi', 'NADP+', 'NADPH', 'THF', '5,10-methylene-THF', '5,10-methenyl-THF', '3-phospho-D-glycerate'};
banRxns = {'human_proteinPool', 'metabolitePool'};

        
[cMatrix, labels, rxnStart, sourceMets, sinkMets] = generateShortestDistance(model, smallSolution, sources, sinks, banMets, banRxns);

graph_to_dot('displayNetwork/test.dot', cMatrix,  labels, rxnStart, sourceMets, sinkMets) 

%Display graph
delete('displayNetwork/test.dot.pdf')
commandStr = 'python displayNetwork/makePdf.py';
 [status, commandOut] = system(commandStr);
 commandOut
 if status==0
     fprintf('pass\n');
 else
     fprintf('fail\n');
 end

end

function graph_to_dot(filename, adj, labels, rxnStart, sourceMets, sinkMets)
    adj = 1000 * adj;

    fid = fopen(filename, 'w');

    fprintf(fid, 'digraph G {\n');
    fprintf(fid, 'rankdir=LR;\n');


    fprintf(fid, 'center = 1;\n');
    fprintf(fid, 'size=\"%d,%d\";\n', 10, 10);


    for i = 1:rxnStart
        additional = ', color = "#E0E0E0", fillcolor = "#E0E0E0", style="filled,setlinewidth(0.1)"';

        if or(sourceMets(i), sinkMets(i))
            additional = ', color = "#5D9BD3", fillcolor = "#5D9BD3", fontcolor = "#FFFFFF", style="filled,setlinewidth(0.1)"';         
        end
        fprintf(fid, '%d [ label = "%s", fontsize=25, shape="box" %s];\n', i, labels{i}, additional);
    end

    allColapsedEdges = [];

    
    
    for i = (rxnStart+1):length(adj)
        letter = char(64+i-rxnStart);
        fprintf(fid, '%d [label = "", shape=circle, fixedsize=true, width=0.1, height=0.1, color = "#A0A0A0", fillcolor = "#A0A0A0", style="filled,setlinewidth(0.1)"];\n', i);
        %fprintf(fid, '%d [label = "%s", shape=circle, fontcolor = "#FFFFFF", fontsize=15, fixedsize=true, width=0.3, height=0.3, fillcolor = "#A0A0A0", style="filled,setlinewidth(0)"];\n', i, letter);
        
    end
    
%Put all sources in a row
sourceIds = find(sourceMets);       
fprintf(fid, '{rank = same; ');
    for i = 1:length(sourceIds)
        fprintf(fid, '%i;', sourceIds(i));
    end
fprintf(fid, '}');
    
    
lineWidth = 2.5;    


%Put all edges

    for node1 = 1:length(adj)
        arcs = find(adj(node1,:));
        
        for node2 = arcs
            flux = adj(node1,node2);
            fprintf(fid, '%d -> %d [arrowsize=1.2, fontsize=8, color="#A0A0A0", penwidth=%f];\n', node1, node2, lineWidth);  
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

