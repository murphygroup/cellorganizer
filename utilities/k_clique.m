function [components,cliques,CC] = k_clique(k,M)
% [X,Y,Z] = k_clique(k,A) 
% Inputs: 
% k - clique size  
% A - adjacency matrix
% 
% Outputs: 
% X - detected communities 
% Y - all cliques (i.e. complete subgraphs that are not parts of larger
% complete subgraphs)
% Z - k-clique matrix

nb_nodes = size(M,1); 
degree_sequence = sort(sum(M,2) - 1,'descend');
max_s = 0;
for i = 1:length(degree_sequence)
    if degree_sequence(i) >= i - 1
        max_s = i;
    else 
        break;
    end
end
cliques = cell(0);
for s = max_s:-1:3
    M_aux = M;
    for n = 1:nb_nodes
        A = n; 
        B = setdiff(find(M_aux(n,:)==1),n); 
        C = transfer_nodes(A,B,s,M_aux); 
        if ~isempty(C)
            for i = size(C,1)
                cliques = [cliques;{C(i,:)}];
            end
        end
        M_aux(n,:) = 0; 
        M_aux(:,n) = 0;
    end
end
CC = zeros(length(cliques));
for c1 = 1:length(cliques)
    for c2 = c1:length(cliques)
        if c1==c2
            CC(c1,c2) = numel(cliques{c1});
        else
            CC(c1,c2) = numel(intersect(cliques{c1},cliques{c2}));
            CC(c2,c1) = CC(c1,c2);
        end
    end
end
CC(eye(size(CC))==1) = CC(eye(size(CC))==1) - k;
CC(eye(size(CC))~=1) = CC(eye(size(CC))~=1) - k + 1;
CC(CC >= 0) = 1;
CC(CC < 0) = 0;
components = [];
for i = 1:length(cliques)
    linked_cliques = find(CC(i,:)==1);
    new_component = [];
    for j = 1:length(linked_cliques)
        new_component = union(new_component,cliques{linked_cliques(j)});
    end
    found = false;
    if ~isempty(new_component)
        for j = 1:length(components)
            if all(ismember(new_component,components{j}))
                found = true;
            end
        end
        if ~found
            components = [components; {new_component}];
        end
    end
end
    function R = transfer_nodes(S1,S2,clique_size,C)
      
        found_s12 = false;
        found_s1 = false;
        for c = 1:length(cliques)
            for cc = 1:size(cliques{c},1)
                if all(ismember(S1,cliques{c}(cc,:)))
                    found_s1 = true;
                end
                if all(ismember(union(S1,S2),cliques{c}(cc,:)))
                    found_s12 = true;
                    break;
                end
            end
        end
        
        if found_s12 || (length(S1) ~= clique_size && isempty(S2))
         
            R = [];
        elseif length(S1) == clique_size
            if found_s1
                R = [];
            else
                R = S1;
            end
        else
            if isempty(find(S2>=max(S1),1))
                R = [];
            else
                R = [];
                for w = find(S2>=max(S1),1):length(S2)
                    S2_aux = S2;
                    S1_aux = S1;
                    S1_aux = [S1_aux S2_aux(w)];
                    S2_aux = setdiff(S2_aux(C(S2(w),S2_aux)==1),S2_aux(w));
                    R = [R;transfer_nodes(S1_aux,S2_aux,clique_size,C)];
                end
            end
        end
    end
end
                
