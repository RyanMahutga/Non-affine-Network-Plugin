function [ fibers_p,deg ] = NetworkPare( fibers, num_nodes, abs_min_deg )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% absolute minimum degree for network

fibers_pared = [];

node1 = fibers(:,1);
node2 = fibers(:,2);

for m = 1:num_nodes
    deg(m) =  sum((node1==m)) + sum((node2==m)) ; % calculate nodal degree
end

init_deg=deg;

list=[];
for m = 1:num_nodes
    list = [find(node1==m); find(node2==m)] ; % indices of fibers connecting to m node
    
    while deg(m) > abs_min_deg && ~isempty(list)
        len = length(list);
        
        r=randi([1,len],1); % randomly select fiber to remove
        
        rem_idx = list(r); % removal index
        if deg(fibers(rem_idx,1)) > abs_min_deg && deg(fibers(rem_idx,2)) > abs_min_deg
            fibers(rem_idx,3) = [NaN]; %remove
            deg(fibers(rem_idx,1)) = deg(fibers(rem_idx,1))-1;
            deg(fibers(rem_idx,2)) = deg(fibers(rem_idx,2)) -1;
        end
        list(r) = [];
    end

end

k=1;
for n = 1:length(fibers)
   if ~isnan(fibers(n,3)) 
       fibers_p(k,:) = fibers(n,:);
       k=k+1;
   end
end

end

