
function [nodes,fibers,fib_rads,init_lens, fib_type, stab_nodes] = periodicNetwork(mu,sigma,scale,num_nodes)
% Create periodic collagen
% Author: Ryan Mahutga
% Barocas Research Group
% Date: 12-03-20

island_scale_factor = 0.95 ; % if this is less than one nodes come from a smaller region (we use this to prevent things being pulled across boundary)

nodes = island_scale_factor.*(rand(num_nodes,3) - 0.5) ;

% select fixed nodes to stabilize networks (closest to [0,0,0])
diststab=[];
% find a node closest to (0,0,0)
for m=1:length(nodes)
    diststab(m) = norm(nodes(m,1:2));
end
[~,indexstab] = min(diststab);
stab_nodes = [indexstab, -1,-1];

num_nodes = length(nodes);

% percolating elastin network
[fibers] = periodicDelaunay(nodes);
[fibers] = removeDupes( fibers );

% [fibers,degC ] = NetworkPare( fibers, length(nodes), 9 );

[ nodes, fibers, stab_nodes ] = MinDegreeReorder( nodes, fibers, stab_nodes ) ; % minimum degree reordering
[ fib_rads, init_lens ] = networkFeatures( nodes, fibers, mu, sigma);

fib_type = 3.*ones(length(fibers));

fiber_vol_fract = sum(pi().*fib_rads.^2.*init_lens*scale)/(scale)^3;
end
