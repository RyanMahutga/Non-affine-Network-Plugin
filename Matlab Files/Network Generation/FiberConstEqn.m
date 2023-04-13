function [ fib_force, node_1_force, node_2_force, dFdlam ] = FiberConstEqn( vect, init_len, fibtype, fiber_area )
%FiberConstEqn contains the fiber constitutive equation for the periodic
%       RVE solution
%   INPUTS: vect - 3x1 vector from node1 to node2
%           init_len - fiber initial length
%           fibtype - fiber type

%   OUTPUTS: fib_force - force of fiber between nodes
%            node_1_force - 3x1 x,y,z force acting on node1
%            node_2_force - 3x1 x,y,z force acting on node2
%            dFdlam - 1x1 value of the derivative of fiber force with
%            respect to lambda at current lambda

%  Created By: Ryan Mahutga - Barocas Research Group - University of Minnesota
%  Date Created: 02/08/18

a = 2; % 1 is for exponential fibers, 2 is for linear fibers, 3 for helical spring model, 4 for incompresible neo-Hookean fibers

% dealing with input values
x_span = vect(1) ;
y_span = vect(2) ;
z_span = vect(3) ;

fiber_current_length = sqrt( x_span^2 + y_span^2 + z_span^2 ); % calculation of fiber length

lambda = fiber_current_length / init_len; % calculation of fiber stretch

% This is an option to use either exponential or linear fibers -- One can
% adjust the fiber parameters below if desired.

%% Exponential Fiber
if a==1  % EXP FIBER FORCE RELATION F = E A ( exp(B*Green Strain) - 1 ) / B
    
    lambda_limit = 7.5; % fib force linear above fiber green strain of 7.5
    
    % generate fiber modulus and beta value
    if fibtype == 1
        fib_mod = 0.7e9; % gives E = 0.5GPa ( %10.8e6; % gives EA=340nN
        fib_B = 2.5;
    elseif fibtype == 2 % elastin fibers
        fib_mod = 2.5e6;
        fib_B = 2.17;
    elseif fibtype == 3 % radial elastin struts (2x modulus of elastin)
        fib_mod = 2.5e6;
        fib_B = 2.17;
    end
    
    if lambda > lambda_limit
        % LINEAR FORCE ABOVE OUR THRESHOLD STRETCH
        
        gs = 0.5 * ( lambda_limit^2 - 1 );
        
        force_exp = fib_mod * fiber_area * (exp(fib_B * gs)-1) / fib_B;
        
        slope_at_lambda_limit = fib_mod * fiber_area * lambda_limit * exp(fib_B * gs);
        
        force = force_exp + slope_at_lambda_limit * (lambda - lambda_limit);
        
        dFdlam = slope_at_lambda_limit ;
        
    elseif lambda < lambda_limit && lambda >= 1
        
        % REGULAR CONST EXN BELOW THRESHOLD
        
        gs = 0.5 * ( lambda^2 - 1 );
        
        force = fib_mod * fiber_area * (exp( fib_B * gs) - 1) / fib_B;
        
        dFdlam = fib_mod * fiber_area * lambda * exp( fib_B * gs) ;
        
    else
        
        %BELOW LAMBDA = 1 -> use 1e-12 * force
        
        gs = 0.5 * ( lambda^2 - 1 );
        
        force = 1e-12.* fib_mod * fiber_area * (exp( fib_B * gs) - 1) / fib_B;
        
        dFdlam = 1e-12* fib_mod * fiber_area * lambda * exp( fib_B * gs) ;
        
%                 force = 0 ;
        
    end
    
    %% Linear Fiber
elseif a==2  % Linear fiber consitutive equation (F=EA*(Green Strain))
    
    % generate fiber modulus
    if fibtype == 1
        fib_mod = 13e6 ; %10.8e6 ; % gives E = 0.5GPa ( %10.8e6; % gives EA=340nN
    elseif fibtype == 2 % elastin fibers
        fib_mod = 1e6;
    elseif fibtype == 3 % radial elastin struts (2x modulus of elastin)
        fib_mod = 170e6;
    end
    
    if lambda >= 1
        
        gs = 0.5 * ( lambda^2 - 1 );
        
        force = fib_mod * fiber_area * gs ;
        
        dFdlam = fib_mod * fiber_area * lambda ;
        
    else % if lambda < 1 (3 orders of magnitude lower
        
        gs = 0.5*(lambda^2-1); % Green Strain
        
        force = 1e-12.*fiber_area*fib_mod*gs ;
        
        dFdlam = 1e-12 * fib_mod * fiber_area * lambda ;
        
    end
    
%% Helical collagen fibers (Freed, A.D. & Doehring, T.C. (2005) J. Biomech. Eng. 127    
elseif a ==3
    
    % generate fiber modulus
    if fibtype == 1
        E = 0.7e9 ; %10.8e6 ; % gives E = 0.5GPa ( %10.8e6; % gives EA=340nN
    elseif fibtype == 2 % elastin fibers
        E = 2.5e6;
    elseif fibtype == 3 % radial elastin struts (2x modulus of elastin)
        E = 2.5e6;
    end
    
    R0 = 2.9 ;% Diameter of 5.8 nm (roughly a collagen microfibril size)
    r0 = 1.2 ; % Diameter of 2.4 nm (roughly a tropocollagen triple helix)
    H0 = 67 ; % [nm] D band pattern in collagen microfibril
    
    L0 = ((2*pi()*R0)^2 + H0^2)^0.5;
    
    lambda_bar = L0/H0 ;
    E_bar = E*H0^2/(L0*(H0 + (1 + 37/(6*pi()^2) + 2*(L0/(pi()*r0))^2)*(L0-H0))) ;
    
    if lambda <= lambda_bar
        H = lambda * H0 ;
        R = (L0^2-H^2)^0.5/(2*pi()) ;
        eta = (R^2+H^2)/(L0*H*(1 + 4*R^2/r0^2 + 6*(20/9 + R^2/r0^2)*R^2/H^2)) ;
        sigma = eta*E_bar*(lambda-1) ;
        dHdlam = H0 ;
        dRdH = -H/(2*pi())*(L0^2-H^2)^(-0.5) ;
        dEtadH = -(H^2+R^2)/(H^2*L0*(6*R^2*(R^2/r0^2+20/9)/H^2 + 4*R^2/r0^2 + 1)) + ...
            2/(L0*(6*R^2*(R^2/r0^2 + 20/9)/H^2 + 4*R^2/r0^2 + 1)) + ...
            12*R^2*(H^2 + R^2)*(R^2/r0^2 + 20/9)/(H^4*L0*(6*R^2*(R^2/r0^2 + 20/9)/H^2 + 4*R^2/r0^2 + 1)^2) ;
        dEtadR = 2*R/(H*L0*(6*R^2*(R^2/r0^2+20/9)/H^2 + 4*R^2/r0^2 + 1)) -...
            (H^2+R^2)*(12*R^2/(r0^2*H^2) + 12*R*(R^2/r0^2 + 20/9)/H^2 + 8*R/r0^2)/(H*L0*(6*R^2*(R^2/r0^2 + 20/9)/H^2 + 4*R^2/r0^2 + 1)^2) ;
        dEtadlam = dEtadH*dHdlam + dEtadR*dRdH*dHdlam ;
        dFdlam = fiber_area*(eta*E_bar + E_bar*(lambda-1)*dEtadlam) ;
    else
        sigma = E_bar*(lambda_bar - 1) + E*(lambda/lambda_bar - 1) ;
        dFdlam = E*fiber_area/lambda_bar ;
    end
    
    force = fiber_area*sigma ; % estimation since sigma is Cauchy Stress
    
%% Neo-Hookean Incompressible Fiber    
elseif a==4 % this is not working currently 
    % Relation obtained assuming 2 and 3 direction are the same, and that 
    % the strain energy density is given by: W = C1/2*(I-3) where I =
    % lam1^2+lam2^2+lam3^2; (Derived in terms of underformed fiber 
    % x-sectional area.)
    
    % generate fiber modulus
    if fibtype == 1
        fib_mod = 30e6 ; %10.8e6 ; % gives E = 0.5GPa ( %10.8e6; % gives EA=340nN
    elseif fibtype == 2 % elastin fibers
        fib_mod = 2.5e6;
    elseif fibtype == 3 % radial elastin struts (2x modulus of elastin)
        fib_mod = 2.5e6;
    end
    
    if lambda >= 1
        
        force = fib_mod * fiber_area * (lambda - 1/(lambda^2)) ;
        
        dFdlam = fib_mod * fiber_area * (1 + 2*(lambda^-3)) ;
        
    else % if lambda < 1 (3 orders of magnitude lower
        
        force = 1e-12 * fib_mod * fiber_area * (lambda - 1/(lambda^2)) ;
        
        dFdlam = 1e-12 * fib_mod * fiber_area * (1 + 2*(lambda^-3)) ;
        
    end
    
end

%% Preparing Outputs
fib_force = force; % fiber forces

cosine_A = x_span / fiber_current_length;
cosine_B = y_span / fiber_current_length;
cosine_C = z_span / fiber_current_length;


% calculating nodal forces

node_1_force(1) = force * cosine_A * (1);
node_1_force(2) = force * cosine_B * (1);
node_1_force(3) = force * cosine_C * (1);

node_2_force(1) = force * cosine_A * (-1);
node_2_force(2) = force * cosine_B * (-1);
node_2_force(3) = force * cosine_C * (-1);

end

