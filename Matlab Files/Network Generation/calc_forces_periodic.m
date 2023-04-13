function [forces, fib_forces, fibers_n, bnd_nodes, nodes_n] = ...
    calc_forces_periodic(nodes, fibers, init_lens, fibtype, fib_areas, rve_size)

%% calc_forces_periodic calculated the forces on nodes and in fibers under periodic boundary conditions

% INPUTS: nodes-Nx3 matrix of internal nodal locations x,y,z
%         fibers-Mx5 matrix of nodes on ends of fibers and fiber crossings
%         [node1 node2 +1 -1 0]
%         init_lens-Mx1 vector of fiber initial lengths
%         modulii-Mx1 vector of fiber moduli (also fib_mods)
%         fiber_area-Mx1 vector of fiber areas
%         fiber_B_vector-Mx1 vecotr of fiber beta values (also fiber_bs)
%         rve_size- 3x3 matrix of rve dimensions[x 0 0; 0 y 0; 0 0 z]

% OUTPUTS: forces-Nx3 matrix of forces on nodes
%          fib_forces-Mx1 matrix of fiber force magnitudes
%          nodes - Nx3 matrix of adjusted internal nodal coordinates (x,y,z)
%          bnd_nodes - Bx3 matrix of boundary node coordinates

% Ryan Mahutga - Barocas Research Group - University of Minnesota
% Date Modified: 03-20-18

%% Calculating fiber forces and nodal forces
num_fibers = length(fibers);
fib_forces = zeros(num_fibers,1);
bnd_nodes = []; % intialize bnd_node matrix
nodes_n = nodes ;

marg = 1e-12; % margin of error built in for bounds

forces = zeros(size(nodes,1),3);

up_bnd(1) = -0.5+rve_size(1,1) + marg ; % x-direction
up_bnd(2) = -0.5+rve_size(2,2) + marg; % y-direction
up_bnd(3) = -0.5+rve_size(3,3) + marg; % z-direction
lo_bnd = -0.5 - marg; % same bound for all directions

for n = 1 : num_fibers

    node1 = fibers(n,1);
    node2 = fibers(n,2);
    
    %% Updating crossing conditions
    for i = 1:3 %(1=x 2=y 3= z)
        current_node1 =  round(nodes(node1,i),12) ;
        current_node2 =  round(nodes(node2,i),12) ;
        k = 1; h = 1; p = 1 ; q = 1;
        while current_node1 > up_bnd(i) || current_node2 > up_bnd(i) || current_node1 < lo_bnd || current_node2 < lo_bnd
            if current_node1 > up_bnd(i)
                fibers(n,2+i) = fibers(n,2+i) - 1;
                nodes_n(node1,i) = nodes(node1,i) - k*rve_size(i,i) ;
                current_node1 = current_node1 - rve_size(i,i) ;
                k = k+1 ;
            elseif current_node1 < lo_bnd
                fibers(n,2+i)= fibers(n,2+i) + 1;
                nodes_n(node1,i) = nodes(node1,i) + p*rve_size(i,i) ;
                current_node1 = current_node1 + rve_size(i,i) ;
                p = p+1;
            end
            
            if current_node2 < lo_bnd
                fibers(n,2+i) = fibers(n,2+i) - 1;
                nodes_n(node2,i) = nodes(node2,i) + h*rve_size(i,i) ;
                current_node2 = current_node2 + rve_size(i,i) ;
                h = h+1;
            elseif current_node2 > up_bnd(i)
                fibers(n,2+i)= fibers(n,2+i) + 1;
                nodes_n(node2,i) = nodes(node2,i) - q*rve_size(i,i) ;
                current_node2 = current_node2 - rve_size(i,i) ;
                q = q+1;
            end
        end
    end
    
    %% Finding Boundary Nodes
    % Generating fiber vectors
    realnode2 = nodes_n(node2,:) + fibers(n,3:5)*rve_size ; % location of node2 on real vector from node1 to node2
    
    x_span = realnode2(1) - nodes_n( node1,1 );
    y_span = realnode2(2) - nodes_n( node1,2 );
    z_span = realnode2(3) - nodes_n( node1,3 );
    
    vect = [x_span, y_span, z_span] ;
    unit = vect./norm(vect) ;
    
    m = 1;
    dl = [];
    for i = 1:3
        if fibers(n,2+i) > 0
            for c = 1:fibers(n,2+i)
                dl(m) = ((-0.5+c*rve_size(i,i)) - nodes_n(node1,i))/unit(i) ; % length of vector to +xi-crossing
                idx(m) = i ;
                m = m+1;
            end
        elseif fibers(n,2+i) < 0
            for c = 1:abs(fibers(n,2+i))
                dl(m) = ((-0.5-(c-1)*rve_size(i,i)) - nodes_n(node1,i))/unit(i) ; % length of vector to -xi-crossing
                idx(m) = i;
                m=m+1;
            end
        end
    end
    
    cross=[];
    [cross, mi] = sort(dl) ;
    if ~isempty(cross)
        for j = 1:length(cross)
            bndnode1 = cross(j).*unit + nodes_n(node1,:) ;
            
            while bndnode1(1) > up_bnd(1) || bndnode1(2) > up_bnd(2) || bndnode1(3) > up_bnd(3)...
                    ||bndnode1(1) < lo_bnd || bndnode1(2) < lo_bnd || bndnode1(3) < lo_bnd
                if bndnode1(1) > up_bnd(1)
                    bndnode1(1) = bndnode1(1) - (rve_size(1,1)) ;
                elseif bndnode1(1) < lo_bnd
                    bndnode1(1) = bndnode1(1) + (rve_size(1,1)) ;
                end
                
                if bndnode1(2) > up_bnd(2)
                    bndnode1(2) = bndnode1(2) - (rve_size(2,2)) ;
                elseif bndnode1(2) < lo_bnd
                    bndnode1(2) = bndnode1(2) + (rve_size(2,2)) ;
                end
                
                if bndnode1(3) > up_bnd(3)
                    bndnode1(3) = bndnode1(3) - (rve_size(3,3)) ;
                elseif bndnode1(3) < lo_bnd
                    bndnode1(3) = bndnode1(3) + (rve_size(3,3)) ;
                end
            end
            
            in = idx(mi(j)) ;
            
            if in == 1 
                A = [1 0 0] ;
            elseif in == 2
                A = [0 1 0] ;
            elseif in == 3
                A = [0 0 1];
            end
            
            if fibers(n,2+in) > 0
                bndnode2 = bndnode1 - A.*rve_size(in,in) ;
            elseif fibers(n,2+in) < 0
                bndnode2 = bndnode1 + A.*rve_size(in,in) ;
            end
            bnd_nodes = [bnd_nodes; bndnode1; bndnode2 ] ;
        end
    end
    
    %% Calculating Nodal and Fiber Forces
    
    [ fib_force, node_1_force, node_2_force, ~ ] = FiberConstEqn( vect, init_lens(n), fibtype(n), fib_areas(n) ) ;
    
    fib_forces(n) = fib_force ; % fiber force
    
    % nodal forces
    forces( node1,1 ) = forces(node1,1) + node_1_force(1);
    forces( node1,2 ) = forces(node1,2) + node_1_force(2);
    forces( node1,3 ) = forces(node1,3) + node_1_force(3);
    
    forces( node2,1 ) = forces(node2,1) + node_2_force(1);
    forces( node2,2 ) = forces(node2,2) + node_2_force(2);
    forces( node2,3 ) = forces(node2,3) + node_2_force(3);
    
    fibers_n = fibers ;
    
end

end
