% Create single type network
% Author: Ryan Mahutga
% Barocas Research Group
% Date: 12-03-20

close all
clear all
clc

plot_net=1;

dir_id='new';
netset=1;

%fiber radii distribution
mu = 30e-9;
sigma = 0e-9;

scale = 10e-6;

num_nodes = 90; % number of elastic lamina nodes

island_scale_factor = 0.9 ; % if this is less than one nodes come from a smaller region (we use this to prevent things being pulled across boundary)

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

 [fibers,degC ] = NetworkPare( fibers, length(nodes), 6);

[ nodes, fibers, stab_nodes ] = MinDegreeReorder( nodes, fibers, stab_nodes ) ; % minimum degree reordering
[ fib_rads, init_lens ] = networkFeatures( nodes, fibers, mu, sigma);

fib_type = ones(length(fibers));

WriteNet3( stab_nodes, fibers, fib_type, init_lens, fib_rads, nodes, scale, dir_id, netset)

fiber_vol_fract = sum(pi().*fib_rads.^2.*init_lens*scale)/(scale)^3;

if plot_net==1
    
    rve_stretch = [1 1 1];
    %fiber_vol_fract = 0.0025;
    
    [ nodes_n, bnd_nodes_n, net_stress, fib_stress, fib_forces, fibers_n ]...
        = solve_periodic_BCs2( nodes, fibers, fib_type, init_lens', fib_rads,...
        fiber_vol_fract, rve_stretch,[] ) ;
    
    %            ConvPeriodic2FixedNetwork( nodes_n, bnd_nodes_n, fibers_n, fib_type)
    
    plot_net_single_fib_type(nodes_n, bnd_nodes_n, fibers_n, fib_type, rve_stretch)
    plot_net_tile_fib_type(nodes_n, bnd_nodes_n, fibers_n, fib_type, rve_stretch, 2)

%     plot3(nodes(stab_nodes(1),1), nodes(stab_nodes(1),2),nodes(stab_nodes(1),3),'bo','LineWidth',3,'MarkerSize',7)
    
end
