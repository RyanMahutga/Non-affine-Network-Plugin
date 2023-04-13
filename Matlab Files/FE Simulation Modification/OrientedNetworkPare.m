function [ fibers_p, deg ] = OrientedNetworkPare(orient, fibers, num_nodes, str )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

abs_min_deg = 6; % absolute minimum degree for network

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
    
    r=0;
    
    hold_orient = orient(list,:);
    if str == "plane"
        [sort_orient,I] = sort(abs(hold_orient(:,3)),'descend');
    elseif str=="dir"
         [sort_orient,I] = sort(abs(hold_orient(:,1)),'ascend');
    else
        disp("ERROR!")
    end
    list2 = list(I);
    len=100;
    while deg(m) > abs_min_deg && ~isempty(list2) && r < len
        len = length(list2);
        r=r+1;
        rem_idx = list2(r); % removal index
        if deg(fibers(rem_idx,1)) > abs_min_deg && deg(fibers(rem_idx,2)) > abs_min_deg
            fibers(rem_idx,3) = [NaN]; %remove
            deg(fibers(rem_idx,1)) = deg(fibers(rem_idx,1))-1;
            deg(fibers(rem_idx,2)) = deg(fibers(rem_idx,2)) -1;
        end
    end
    list2=[];
    list=[];
    hold_orient=[];
    sort_orient=[];
end

% node1=[];
% node2=[];
% 
% node1 = fibers(:,1);
% node2 = fibers(:,2);
% 
% % recalculate degree
% for m = 1:num_nodes
%     deg(m) =  sum((node1==m)) + sum((node2==m)) ; % calculate nodal degree
% end

k=1;
for n = 1:length(fibers)
   if ~isnan(fibers(n,3)) 
       fibers_p(k,:) = fibers(n,:);
       k=k+1;
   end
end

end

