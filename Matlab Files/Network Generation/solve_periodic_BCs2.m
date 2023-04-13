function [ nodes_n, bnd_nodes_n, net_stress, fib_stress, fib_forces, fibers_n ]...
    = solve_periodic_BCs2( nodes,fibers, fibtype, init_lens, fib_rads, ...
    fiber_vol_fract, rve_stretch, rel_stretch )
%solve_periodic solves the network problem having periodic boundary
%   conditions. NOTE: Network must be periodic for this solution to be
%   applied appropriately.

%   INPUTS: nodes-Nx3 matrix of internal nodal coordinates
%           fibers-Mx5 matrix of nodal connectivity of fibers and crossing
%           boundaries [node1 node2 +1 -2 0] (+1 indicates crosses positive
%           boundary, -1 indicates crosses negative boundary twice, i.e. the
%           fiber wraps around the rve, 0  indicates no crossing of boundary)
%           fibtype-Mx1 vector of fiber types
%           init_lens-Mx1 vector of fiber initial lengths
%           fib_rads-Mx1 vector of fiber radii
%           fiber_vol_fract-1x1 value of fiber volume fraction
%           rve_stretch-3x1 vector of network stretch in x,y,z

% OUTPUTS: nodes_n-Nx3 matrix of deformed internal nodal coordinates (x,y,z)
%          bnd_nodes_n-Bx3 matrix of boundary nodes after solution
%          net_stress-3x1 vector of stress computed from calc_net_stress
%          fib_stress-Nx1 vector of fiber stress
%          fib_forces-Nx1 vector of fiber forces

% Ryan Mahutga - Barocas Research Group - University of Minnesota
% Date Modified: 02-08-18

%% Dealing with inputs and setting up parameters
num_nodes = size(nodes,1) ;
num_fibers = size(fibers,1) ;

fib_areas = pi().*fib_rads.^2 ;

% stretch parameters
x_stretch = rve_stretch(1);
y_stretch = rve_stretch(2);
z_stretch = rve_stretch(3);

rve_size = [x_stretch, 0, 0 ;
            0, y_stretch, 0 ;
            0, 0, z_stretch] ;

if isempty(rel_stretch) % stretch net box -- assumes RVE origin corner is fixed at (-0.5,-0.5,-0.5)
    nodes(:,1) = ((nodes(:,1)+0.5) * x_stretch) - 0.5;
    nodes(:,2) = ((nodes(:,2)+0.5) * y_stretch) - 0.5;
    nodes(:,3) = ((nodes(:,3)+0.5) * z_stretch) - 0.5;
else
    nodes(:,1) = ((nodes(:,1)+0.5) * rel_stretch(1)) - 0.5;
    nodes(:,2) = ((nodes(:,2)+0.5) * rel_stretch(2)) - 0.5;
    nodes(:,3) = ((nodes(:,3)+0.5) * rel_stretch(3)) - 0.5;
end

% extracting free node coordinates
X0 = reshape(nodes', 1, []);

%% Solving for new internal nodal positions (Newton's Method)
for n = 1 : 251 % GIVE UP AT 250 ITERS -- increase this to get better convergence (10^-12 is pretty good)
    
    nodes_X0(:,1) = X0(1:3:3*num_nodes) ;
    nodes_X0(:,2) = X0(2:3:3*num_nodes) ;
    nodes_X0(:,3) = X0(3:3:3*num_nodes) ;
    
    nodes0 = nodes_X0;
    
    [forces, ~, ~, ~, ~] = calc_forces_periodic(nodes_X0, fibers, init_lens, fibtype, fib_areas, rve_size) ;
    
    residuals = reshape(forces', 1, []);
    
    err = norm(residuals);
    
    if err < 1e-11 % fixed newton loop tolerance
        % or INITIAL_ERR_VAL * 1e-5 ==> relative error
        %         disp('Converged')
        break
    end
    
    [ Jac ] = calc_Jacobian2( nodes_X0, fibers, init_lens, fibtype, fib_areas, rve_size  ) ;
    
    %disp(cond(Jac))
    
    J = sparse(Jac) ;
    
    % using the gmres solver
    
    if 1 == 1
        % [L,U] = ilu(sparse(J),struct('type','ilutp','droptol',1e-12));
        [delta,~,~] = gmres(J, -residuals', 10, 1e-12, 100); % 50-30 works as well
    end
    
    % using the backslash solver
    if 1 == 0
        delta = J \ -residuals'; % transpose fvec to column vector
        if n < 50
            delta = 0.10 * delta;
        end
    end
    
    if n==251
        disp(strcat('Failed to Converge: Error=', num2str(err)))
    end
    
    %disp(cond(Jac))
    
    X0 = X0 + delta';
   
end

%% Generating final results

% apply X0 to nodal coordinates
nodes(:,1) = X0(1:3:3*num_nodes) ;
nodes(:,2) = X0(2:3:3*num_nodes) ;
nodes(:,3) = X0(3:3:3*num_nodes) ;

[forces, fib_forces, fibers_n, bnd_nodes_n, nodes_n] = calc_forces_periodic(nodes, fibers, init_lens, fibtype, fib_areas, rve_size) ;
fib_stress = fib_forces./fib_areas ;

[stress] = calc_net_stress_periodic(nodes_n, bnd_nodes_n, fibers_n, rve_size, fib_forces) ;

real_fib_vol = sum(init_lens.*fib_areas') ;
real_stress_factor = ( fiber_vol_fract / real_fib_vol ) ;

current_x_length = rve_size(1,1);
current_y_length = rve_size(2,2);
current_z_length = rve_size(3,3); %when collapsing, this will be a negative number so it will subtract

new_volume = current_x_length * current_y_length * current_z_length;
 
vol_term = 1 / new_volume;

net_stress = vol_term*real_stress_factor.*stress ;

end

