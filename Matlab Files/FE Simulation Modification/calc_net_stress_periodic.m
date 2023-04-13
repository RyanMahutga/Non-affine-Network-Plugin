function [stress] = calc_net_stress_periodic(nodes, bnd_nodes, fibers, rve_size, fib_forces)

% calc_net_stress_periodic.m calculates the volume averaged network stress
%       for a periodic network

% INPUTS: nodes - Nx3 array of interior node coordinates (xyz)
%         bnd_node_nums - Bx3 array of boundary node coordinates (xyz)
%         fibers - Mx5 matrix of nodal coordinates (columns 1 and 2) and 
%                   crossing data (x,y,z in columns 3,4,5) 
%         fib_forces - Mx1 array of fiber forces between nodes
%
% OUTPUTS: stress -- 3x3 array for the volume averaged stress tensor

% Created by: Ryan Mahutga - Barocas Research Group - University of Minnesota 
% Date Modified: 02-01-18

% Finding Fiber Crossing Forces

bnd_coords = zeros(length(bnd_nodes),3) ;
bnd_forces = zeros(length(bnd_nodes),3) ;
m = 1;
for n = 1:length(fibers)
    node1 = fibers(n,1) ;
    node2 = fibers(n,2) ;
    
    % Generating fiber vectors
    realnode2 = nodes(node2,:) + fibers(n,3:5)*rve_size ; % location of node2 on real vector from node1 to node2
    
    x_span = realnode2(1) - nodes( node1,1 );
    y_span = realnode2(2) - nodes( node1,2 );
    z_span = realnode2(3) - nodes( node1,3 );
    
    vect = [x_span, y_span, z_span] ;
    fiber_current_length = sqrt( x_span^2 + y_span^2 + z_span^2 );  
    unit = vect./fiber_current_length ;
    
       for i = 1:3
            for k = 1:abs(fibers(n,2+i))
                bnd_coords(m,:) = bnd_nodes(m,:) ;                
                bnd_forces(m,:) = unit.*fib_forces(n) ;
                m = m+1 ;
                
                bnd_coords(m,:) = bnd_nodes(m,:) ;
                bnd_forces(m,:) = -unit.*fib_forces(n) ;
                m = m+1 ;
            end
       end

end

xmax = max(bnd_coords(:,1));
ymax = max(bnd_coords(:,2));
zmax = max(bnd_coords(:,3));
mini = -0.5 ;

xfx=0; yfx=0; zfx=0;
xfy=0; yfy=0; zfy=0;
xfz=0; yfz=0; zfz=0;

% Sum xi * fi through matrix multiplication
xfx = (bnd_coords(:,1)' * bnd_forces(:,1));
yfx = (bnd_coords(:,2)' * bnd_forces(:,1));
zfx = (bnd_coords(:,3)' * bnd_forces(:,1));

xfy = (bnd_coords(:,1)' * bnd_forces(:,2));
yfy = (bnd_coords(:,2)' * bnd_forces(:,2));
zfy = (bnd_coords(:,3)' * bnd_forces(:,2));

xfz = (bnd_coords(:,1)' * bnd_forces(:,3));
yfz = (bnd_coords(:,2)' * bnd_forces(:,3));
zfz = (bnd_coords(:,3)' * bnd_forces(:,3));

stress = [xfx, yfx, zfx; xfy, yfy, zfy; xfz, yfz, zfz];

% Formula for this is Sij=1/V*Sum(xiFj); The off diagnonal terms are
% symmetric, e.g. xfy=yfx. However, there is a very small difference in
% these numbers due to numerical issues, so we just average them. The
% difference seems to be at the 10-14 decimal place, so probably could
% ignore.

% compact notation [S11 S12 S13 S22 S23 S33]
%stress= [xfx; 0.5*(xfy+yfx); 0.5*(zfx+xfz); yfy;0.5*(yfz+zfy); zfz]; 

% full notation [S11 S12 S13; S21 S22 S23; S31 S32 S33]
% stress= [ xfx 0.5*(xfy+yfx) 0.5*(zfx+xfz); 
%           0.5*(xfy+yfx) yfy 0.5*(yfz+zfy); 
%           0.5*(zfx+xfz) 0.5*(yfz+zfy) zfz ];

end
