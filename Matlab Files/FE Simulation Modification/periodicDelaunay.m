function [fibers] = periodicDelaunay(nodes)
% periodicDelaunay creates a delaunay tesselation for a periodic network
%   Detailed explanation goes here

    num_nodes = length(nodes);

 % vector to shift nodal coordinates
    add_vect = ...
        [0 0 1; 0 0 -1;
        0 1 0; 0 -1 0;
        1 0 0; -1 0 0;
        1 1 0; 1 -1 0; -1 1 0; -1 -1 0;
        0 1 1; 0 -1 1; 0 1 -1; 0 -1 -1;
        1 0 1; -1 0 1; 1 0 -1; -1 0 -1;
        1 1 1; -1 1 1; 1 -1 1; 1 1 -1;
        -1 -1 1; -1 1 -1; 1 -1 -1; -1 -1 -1] ;
    
    % construct correspondence vector
    correspondence = [];
    for n = 1:26
        for j=1:num_nodes
            correspondence = [correspondence ; j, add_vect(n,:)];
        end
        % generate nodes outside RVE in 3x3 space
        out_nodes((n-1)*num_nodes+1:n*num_nodes ,:) =...
            nodes + add_vect(n,:) ;
    end
    
    all_nodes = [nodes;out_nodes]; % all nodes

    % generating fibers (Delaunay)
    TRI = delaunay(all_nodes(:,1),all_nodes(:,2),all_nodes(:,3));
    
    [row,col] = find(TRI<=num_nodes);
    fibers=[];
    for j=1:length(row)
        if col(j) == 4
            if TRI(row(j),col(j)-1) > num_nodes
                idx = TRI(row(j),col(j)-1) - num_nodes ;
                fibers = [fibers; TRI(row(j),col(j)), correspondence(idx,:)];
            else
                fibers = [fibers; TRI(row(j),col(j)), TRI(row(j),col(j)-1),0,0,0];
            end
        elseif col(j) == 1
            if TRI(row(j),col(j)+1) > num_nodes
                idx = TRI(row(j),col(j)+1) - num_nodes ;
                fibers = [fibers; TRI(row(j),col(j)), correspondence(idx,:)];
            else
                fibers = [fibers; TRI(row(j),col(j)), TRI(row(j),col(j)+1),0,0,0];
            end
        else
            if TRI(row(j),col(j)+1) > num_nodes
                idx = TRI(row(j),col(j)+1) - num_nodes ;
                fibers = [fibers; TRI(row(j),col(j)), correspondence(idx,:)];
            else
                fibers = [fibers; TRI(row(j),col(j)), TRI(row(j),col(j)+1),0,0,0];
            end
            
            if TRI(row(j),col(j)-1) > num_nodes
                idx = TRI(row(j),col(j)-1) - num_nodes ;
                fibers = [fibers; TRI(row(j),col(j)), correspondence(idx,:)];
            else
                fibers = [fibers; TRI(row(j),col(j)), TRI(row(j),col(j)-1),0,0,0];
            end
        end
    end
end

