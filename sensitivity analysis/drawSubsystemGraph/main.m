labelsA = {'Nucleotide metabolism'
'Valine, leucine, and isoleucine metabolism'
'Glycolysis / Gluconeogenesis'
'Fatty acid biosynthesis and transfer reactions'
'Oxidative phosphorylation'
'Cholesterol biosynthesis and metabolism'
'Cysteine and methionine metabolism'
'Starch and sucrose metabolism'
'Alanine, aspartate and glutamate metabolism'
};
labelsB = {'Membrane'; 'Nucleotides'; 'Protein'; 'Energy'};
labelsC = {'Growth'};
labels= [labelsA;labelsB;labelsC];

adj = zeros(length(labels));
adj(findIndex(labels, 'Nucleotide metabolism'), findIndex(labels, 'Nucleotides')) = 1;
adj(findIndex(labels, 'Valine, leucine, and isoleucine metabolism'), findIndex(labels, 'Energy')) = 1;
adj(findIndex(labels, 'Valine, leucine, and isoleucine metabolism'), findIndex(labels, 'Membrane')) = 1;
adj(findIndex(labels, 'Glycolysis / Gluconeogenesis'), findIndex(labels, 'Energy')) = 1;
adj(findIndex(labels, 'Fatty acid biosynthesis and transfer reactions'), findIndex(labels, 'Membrane')) = 1;
adj(findIndex(labels, 'Oxidative phosphorylation'), findIndex(labels, 'Energy')) = 1;
adj(findIndex(labels, 'Cholesterol biosynthesis and metabolism'), findIndex(labels, 'Membrane')) = 1;
%adj(findIndex(labels, 'Cysteine and methionine metabolism'), findIndex(labels, 'Membrane')) = 1;
adj(findIndex(labels, 'Starch and sucrose metabolism'), findIndex(labels, 'Energy')) = 1;
adj(findIndex(labels, 'Alanine, aspartate and glutamate metabolism'), findIndex(labels, 'Nucleotides')) = 1;
adj(findIndex(labels, 'Alanine, aspartate and glutamate metabolism'), findIndex(labels, 'Protein')) = 1;


adj(findIndex(labels, 'Membrane'), findIndex(labels, 'Growth')) = 1;
adj(findIndex(labels, 'Nucleotides'), findIndex(labels, 'Growth')) = 1;
adj(findIndex(labels, 'Protein'), findIndex(labels, 'Growth')) = 1;
adj(findIndex(labels, 'Energy'), findIndex(labels, 'Growth')) = 1;







fid = fopen('test.dot', 'w');
fprintf(fid, 'digraph G {\n');
fprintf(fid, 'center = 1;\n');
fprintf(fid, 'rankdir=LR;\n');
fprintf(fid, 'rank = same;\n');

fprintf(fid, 'graph [pad="0", nodesep="0.1", ranksep="2"];\n');

fprintf(fid, 'size=\"%d,%d\";\n', 10, 10);

A=length(labelsA);
for i = 1:A
    fprintf(fid, '%d [ label = "%s\\r", pos="1,%2.2f!", fontsize=20, height=0.3, width=5.4, shape="box", fillcolor = "#97B9E0", style="filled,setlinewidth(0)"];\n', i, labels{i},  i);
end

B = length(labelsB);
for i = (A+1):(A+B)
    fprintf(fid, '%d [ label = "%s\\l", fontsize=20, height=0.3, width=1.6, shape="box", fillcolor = "#97B9E0", style="filled,setlinewidth(0)"];\n', i, labels{i});
end

C = length(labelsC);
for i = (A+B+1):(A+B+C)
    fprintf(fid, '%d [ label = "%s\\r", fontsize=20, height=0.3, shape="box", fillcolor = "#97B9E0", style="filled,setlinewidth(0)"];\n', i, labels{i});
end


lineWidth = 1.8;    



for node1 = 1:length(adj)
    arcs = find(adj(node1,:));

    for node2 = arcs
        flux = adj(node1,node2);
        fprintf(fid, '%d -> %d [arrowsize=0.5, color="#A0A0A0", penwidth=%f];\n', node1, node2, lineWidth);  
    end
end    
fprintf(fid, '}\n'); 
fclose(fid);
