function [omega, theta, phi] = calcOrientation(nodes, fibers, fib_rads, F)
%caclulates the oreitnation tensor omega, and the individual fiber
%orientations theta and phi.

omega = zeros(3,3);

tot_fib_vol = 0.0;

for f=1:size(fibers,1)
    n1 = fibers(f,1);
    n2 = fibers(f,2);

    real_node2(1) = nodes(n2,1) + F(1,1)*fibers(f,3) + F(1,2)*fibers(f,4) + F(1,3)*fibers(f,5) ;
    real_node2(2) = nodes(n2,2) + F(2,1)*fibers(f,3) + F(2,2)*fibers(f,4) + F(2,3)*fibers(f,5) ;
    real_node2(3) = nodes(n2,3) + F(3,1)*fibers(f,3) + F(3,2)*fibers(f,4) + F(3,3)*fibers(f,5) ;

    node1 = nodes(n1,:);

    L = real_node2 - node1;

    % calculating orientation tensor omega (normalized by total fiber
    % volume)
    cosa = L(1)/norm(L);
    cosb = L(2)/norm(L);
    cosg = L(3)/norm(L);

    fib_vol = fib_rads(f)*norm(L);

    omega(1,1) = omega(1,1) + fib_vol * cosa * cosa;
    omega(1,2) = omega(1,2) + fib_vol * cosa * cosb;
    omega(1,3) = omega(1,3) + fib_vol * cosa * cosg;
    omega(2,2) = omega(2,2) + fib_vol * cosb * cosb;
    omega(2,3) = omega(2,3) + fib_vol * cosb * cosg;
    omega(3,3) = omega(3,3) + fib_vol * cosg * cosg;

    tot_fib_vol = tot_fib_vol + fib_vol;

    % calculating single theta/phi
    theta(f) = atand(L(2)/L(1));
    phi(f) = atand(L(3)/(sqrt(L(1)^2+L(2)^2)));
end

% fill in the symmetric values
omega(2,1) = omega(1,2);
omega(3,1) = omega(1,3);
omega(3,2) = omega(2,3);

omega = omega./tot_fib_vol; % normalize by volume

end