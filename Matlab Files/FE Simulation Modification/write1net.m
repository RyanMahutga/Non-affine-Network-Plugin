
close all
clear
clc

mu=150e-9;
sigma=10e-9;

scale = 10e-6;
num_nodes = 25;

dir_id = "singleNetCollagen";
mkdir(dir_id)
netset=1;

[nodes,fibers,fib_rads,init_lens, fib_type, stab_nodes] = periodicNetwork(mu,sigma,scale,num_nodes);

fiber_vol_fract = sum(pi().*fib_rads.^2.*init_lens)/(scale^2);

WriteNet3( stab_nodes, fibers, fib_type, init_lens, fib_rads, nodes, scale, dir_id, netset)

rve_stretch = [1 1 1];

[ nodes_n, bnd_nodes_n, net_stress, fib_stress, fib_forces, fibers_n ]...
    = solve_periodic_BCs2( nodes, fibers, fib_type, init_lens', fib_rads,...
    fiber_vol_fract, rve_stretch,[] ) ;

%            ConvPeriodic2FixedNetwork( nodes_n, bnd_nodes_n, fibers_n, fib_type)

plot_net_single_fib_type(nodes_n, bnd_nodes_n, fibers_n, fib_type, rve_stretch)

plot3(nodes(stab_nodes(1),1), nodes(stab_nodes(1),2),nodes(stab_nodes(1),3),'bo','LineWidth',3,'MarkerSize',7)

