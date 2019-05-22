function linkAA(model, smallSolution, sources, sinks, banMets, banRxns)

for i = 1:length(model.rxns)
    subRxns = split(model.rxns{i},'+');
    model.rxns{i} = subRxns{1};
end

essential = {'histidine'
            'isoleucine'
            'leucine'
            'lysine'
            'methionine'
            'phenylalanine'
            'threonine'
            'tryptophan'
            'valine'};
        
[cMatrix, labels, rxnStart, sourceMets, sinkMets] = generateShortestDistance(model, smallSolution, sources, sinks, banMets, banRxns);

fileName = 'displayNetwork/test.dot';

graph_to_dot(fileName, cMatrix,  labels, rxnStart, sourceMets, sinkMets, essential) 

%Display graph
fileEnding = 'pdf'; %pdf or eps
delete(['displayNetwork/test.dot.' fileEnding])
fileName = [cd '\' fileName];
fileName = strrep(fileName, '\', '/');

commandStr = ['python displayNetwork/makePdf.py ' fileName ' ' fileEnding];
 [status, commandOut] = system(commandStr);
 
 if status==0
     fprintf('pass\n');
 else
     fprintf('fail\n');
 end

end

function graph_to_dot(filename, adj, labels, rxnStart, sourceMets, sinkMets, essential)
    adj = 1000 * adj;

    fid = fopen(filename, 'w');

    fprintf(fid, 'digraph G {\n');
    fprintf(fid, 'rankdir=LR;\n');


    fprintf(fid, 'center = 1;\n');
    fprintf(fid, 'size=\"%d,%d\";\n', 10, 10);


    for i = 1:rxnStart
        %base node:
        additional = ', color = "#E0E0E0", fillcolor = "#E0E0E0"';
        
        if sourceMets(i)
            additional = ', color = "#84C557", fillcolor = "#84C557", fontcolor = "#FFFFFF"';
        elseif sinkMets(i)
            additional = ', color = "#9876A8", fillcolor = "#9876A8", fontcolor = "#FFFFFF"';    
        end      

        additional = [additional ', style="filled,setlinewidth(0.1)"'];
        fprintf(fid, '%d [ label = "%s", fontsize=25, shape="box" %s];\n', i, labels{i}, additional);
    end

    
    
    allColapsedEdges = [];

    for i = (rxnStart+1):length(adj)
        letter = char(64+i-rxnStart);
        if onlyOneEdge(adj, i)
            allColapsedEdges = [allColapsedEdges i];
        else
            fprintf(fid, '%d [label = "", shape=circle, fixedsize=true, width=0.1, height=0.1, color = "#A0A0A0", fillcolor = "#A0A0A0", style="filled,setlinewidth(0.1)"];\n', i);
            %fprintf(fid, '%d [label = "%s", shape=circle, fontcolor = "#FFFFFF", fontsize=15, fixedsize=true, width=0.3, height=0.3, fillcolor = "#A0A0A0", style="filled,setlinewidth(0)"];\n', i, letter);
        end
    end
    
    
%Put all sources in a row
sourceIds = find(sourceMets);       
fprintf(fid, '{rank = same; ');
    for i = 1:length(sourceIds)
        fprintf(fid, '%i;', sourceIds(i));
    end
fprintf(fid, '}');
    
    
lineWidth = 1.5;    


%Draw all edges

    for node1 = 1:length(adj)
        arcs = find(adj(node1,:));
        
        if not(ismember(node1, allColapsedEdges))
            for node2 = arcs
                if ismember(node2, allColapsedEdges)
                    node2 = find(adj(node2,:));                    
                end
                flux = adj(node1,node2);
                fprintf(fid, '%d -> %d [arrowsize=1.5, fontsize=8, color="#A0A0A0", penwidth=%f];\n', node1, node2, lineWidth);  
            end
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

function pass = onlyOneEdge(adj, i)
    connectedNodes1 = sum(adj(i,:)~=0);
    connectedNodes2 = sum(adj(:,i)~=0);
    pass = and(connectedNodes1 == 1, connectedNodes2 == 1);
end

