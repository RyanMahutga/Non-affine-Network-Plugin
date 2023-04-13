function [ J ] = calc_Jacobian2( nodes, fibers, init_lens, fibtype, fib_areas, rve_size  )
% calc_Jacobian calculates the numerical Jacobian of the nodal forces
%   w.r.t. perturbations to the nodal x,y,z coordinates. The Jacobian
%   is calculated as: J = d(F^i)/d(x^j) = d(mag(F)*n^i)/d(x^j) =
%   d(mag(F^i))/d(lambda^i) * d(lambda^i)/dx^j * n^i + d(n^i)/d(x^j) * mag(F^i).
%   F is the force between two fiber nodes, n^i is the normal vector
%   contribution in a given direction between two fiber, x^j is the
%   coordinate being probed, lambda is the fiber stretch between two nodes.
%   All derivative quantities are calculated numerically for dx prescribed
%   below.

% NOTE: The value of dx (aka epsilon) greatly affects the calculation of
%       the numerical Jacobian!

%   INPUTS: nodes - Nx3 matrix of nodal x,y,z coordinates
%           fibers - Mx5 vector of fiber connectivity
%           init_lens - Mx1 vector of fiber initial lengths
%           fibtype - Mx1 vector of fiber types
%           fib_areas - Mx1 vector of fiber cross-sectional areas
%           rve_size - 3x1 current rve size given the stretch of the
%           network

%   OUTPUTS: J - 3M x 3M Jacobian matrix dF^i/dx^j

% defining dx (epsilon)
dx = 1e-8 ;

J = zeros(3*size(nodes,1)) ;
nodes_n = nodes;

marg = 1e-12; % margin of error for bounds
up_bnd(1) = -0.5+rve_size(1,1)+marg ; % x-direction
up_bnd(2) = -0.5+rve_size(2,2)+marg ; % y-direction
up_bnd(3) = -0.5+rve_size(3,3)+marg ; % z-direction
lo_bnd = -0.5-marg ; % same bound for all directions

% tether force Jacobian (added for stability)
% NOTE: There is no dFi/dxj since we assume the force on one node is
% unaffected by the position of another.

% k = 1e-12 ; % proportionality constant for tether (arbitrary)
% for m = 1:length(nodes)
%     % dFi/dxi
%     J(3*m-2 , 3*m-2) = -2*k ;
%     J(3*m-1 , 3*m-1) = -2*k ;
%     J(3*m-0 , 3*m-0) = -2*k ;
%     
% end

for n = 1:length(fibers)
    
    node1 = fibers(n,1) ;
    node2 = fibers(n,2) ;
    
    % shifting node2 to real position
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
    
    realnode2 = nodes_n(node2,:) + fibers(n,3:5)*rve_size ; % location of node2 on real vector from node1 to node2
    
    % calculating original values
    vect = realnode2-nodes_n(node1,:) ;
    unit = vect./norm(vect) ; % fiber unit vector
    
    len = norm(vect);
    
    lambda = norm(vect)/init_lens(n) ;
    
    [ ~, node_1_force, node_2_force, dFdlam] = FiberConstEqn( vect, init_lens(n), fibtype(n), fib_areas(n) ) ;
           
    % node 1 and 2 force should be the same (opposite directions), we average them to mitigate
    % any slight differences
    force = (norm(node_1_force) + norm(node_2_force))/2;
    
     dlamdx1 = -1/(init_lens(n)*len).*vect;
     dlamdx2 = 1/(init_lens(n)*len).*vect;
     
     for i=1:3
         for j=1:3
             
             if i==j
                 dndx11 = -1/len + vect(i)^2/(len^3) ;
                 dndx12 = 1/len - vect(i)^2/(len^3) ;
                 dndx21 = -1/len + vect(i)^2/(len^3) ;
                 dndx22 = 1/len - vect(i)^2/(len^3) ;
             else
                 dndx11 = vect(i)*vect(j)/(len^3) ;
                 dndx12 = -vect(i)*vect(j)/(len^3) ;
                 dndx21 = vect(i)*vect(j)/(len^3) ;
                 dndx22 = -vect(i)*vect(j)/(len^3) ;
             end
                          
             % Calculate dFi/dxj
             j_temp11 = dFdlam * dlamdx1(j) * unit(i) + force * dndx11 ; 
             j_temp12 = dFdlam * dlamdx2(j) * unit(i) + force * dndx12 ; 
             j_temp21 = -dFdlam * dlamdx1(j) * unit(i) - force * dndx21 ; % unit vector (n) has opposite sign
             j_temp22 = -dFdlam * dlamdx2(j) * unit(i) - force * dndx22 ; % unit vector (n) has opposite sign
             
             % this is dFi/dxi
             J(3*node1-(3-i) , 3*node1-(3-j)) = J(3*node1-(3-i) , 3*node1-(3-j)) + j_temp11 ;
             % this is dFj/dxj
             J(3*node2-(3-i) , 3*node2-(3-j)) = J(3*node2-(3-i) , 3*node2-(3-j)) + j_temp22 ;
             
             % this is dFi/dxj
             J(3*node1-(3-i) , 3*node2-(3-j)) = J(3*node1-(3-i) , 3*node2-(3-j)) + j_temp12 ;
             % this is dFj/dxi
             J(3*node2-(3-i) , 3*node1-(3-j)) = J(3*node2-(3-i) , 3*node1-(3-j)) + j_temp21 ;
         end
     end
     
    
end

end

