function plot_net_tile_fib_type(nodes, bnd_nodes, fibers, fib_type, rve_stretch, plane)

% plot_net_tile.m plots a 3x3 repeated RVE in the plane specifed (x=1, y=2,z=3)
%   INPUTS: nodes-Nx3 matrix of deformed internal nodal coordinates (x,y,z)
%           bnd_nodes-Bx3 matrix of boundary nodes after solution
%           fibers-Mx5 matrix of nodal connectivity of fibers and crossing
%           boundaries [node1 node2 +1 -2 0] (+1 indicates crosses positive
%           boundary, -1 indicates crosses negative boundary twice, i.e. the
%           fiber wraps around the rve, 0  indicates no crossing of boundary)
%           rve_stretch-3x1 vector of network stretch in x,y,z
%           plane - 1x1 value indicating plane of viewing (1=x, 2=y, 3=z)
%   OUTPUTS: No variable outputs. Figure should display the 3x3 repeated RVE

% Created by: Ryan Mahutga - Barocas Research Group - University of Minnesota
% Date created: 02-15-18

% check plane input
if plane == 1
    g(1) = 2 ;
    g(2) = 3 ;
elseif plane == 2
    g(1) = 1 ;
    g(2) = 3 ;
elseif plane == 3
    g(1) = 1 ;
    g(2) = 2 ;
else
    disp('Plane must be 1 for x, 2 for y, and 3 for z!')
    return
end

%% Plotting central RVE
figure;

plot3(nodes(:,1) , nodes(:,2), nodes(:,3),'o', 'LineWidth',0.2, 'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',5)
hold on
plot3(bnd_nodes(:,1) , bnd_nodes(:,2), bnd_nodes(:,3),'X', 'LineWidth',0.4, 'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',5);

m = 1 ;

for n = 1 : length(fibers) % count rows
    
    node1 = nodes(fibers(n,1),:) ;
    
    for i = 1:3
        for k = 1:abs(fibers(n,2+i))
            node2 = bnd_nodes(m,:) ;
            
            x1(1) = node1(1); % node 1 x coord
            y1(1) = node1(2); % node 1 y coord
            z1(1) = node1(3); % node 1 z coord
            
            x1(2) = node2(1); % node 2 x coord
            y1(2) = node2(2); % node 2 y coord
            z1(2) = node2(3); % node 2 z coord
            
            if fib_type(n) == 1
                C = [0.8 0.8 0.4];
                lw = 1.5;
            elseif fib_type(n) == 2
                C = [0 0.2 0];
                lw=2;
            elseif fib_type(n) == 3
                C = [0.7 0 0];
                lw=1.5;
            end
            
            plot3(x1, y1, z1, 'Color', C,'LineWidth',lw);
            
            m = m+1 ;
            
            node1 = bnd_nodes(m,:) ;
            
            m = m+1 ;
        end
    end
    
    node2 = nodes(fibers(n,2),:) ;
    
    x1(1) = node1(1); % node 1 x coord
    y1(1) = node1(2); % node 1 y coord
    z1(1) = node1(3); % node 1 z coord
    
    x1(2) = node2(1); % node 2 x coord
    y1(2) = node2(2); % node 2 y coord
    z1(2) = node2(3); % node 2 z coord
    
    if fib_type(n) == 1
        C = [0.8 0.8 0.4];
        lw = 1.5;
    elseif fib_type(n) == 2
        C = [0 0.2 0];
        lw = 2;
    elseif fib_type(n) == 3
        C = [0.7 0 0];
        lw = 1.5;
    end
    
    plot3(x1, y1, z1, 'Color', C, 'LineWidth',lw);
end

%% Plotting the 8 surrounding RVEs
vect1 = ones(length(nodes),1) ;
vect2 = ones(length(bnd_nodes),1);

for p = 1:8
    nodesf = nodes ;
    bnd_nodesf = bnd_nodes ;
    if p == 1
        nodesf(:,g(1)) = nodes(:,g(1)) + vect1*rve_stretch(g(1)) ;
        bnd_nodesf(:,g(1)) = bnd_nodes(:,g(1)) + vect2*rve_stretch(g(1)) ;
    elseif p == 2
        nodesf(:,g(1)) = nodes(:,g(1)) + vect1*rve_stretch(g(1)) ;
        nodesf(:,g(2)) = nodes(:,g(2)) + vect1*rve_stretch(g(2)) ;
        bnd_nodesf(:,g(1)) = bnd_nodes(:,g(1)) + vect2*rve_stretch(g(1)) ;
        bnd_nodesf(:,g(2)) = bnd_nodes(:,g(2)) + vect2*rve_stretch(g(2)) ;
    elseif p == 3
        nodesf(:,g(2)) =  nodes(:,g(2)) + vect1*rve_stretch(g(2)) ;
        bnd_nodesf(:,g(2)) = bnd_nodes(:,g(2)) + vect2*rve_stretch(g(2)) ;
    elseif p == 4
        nodesf(:,g(1)) = nodes(:,g(1)) - vect1*rve_stretch(g(1)) ;
        nodesf(:,g(2)) = nodes(:,g(2)) + vect1*rve_stretch(g(2)) ;
        bnd_nodesf(:,g(1)) = bnd_nodes(:,g(1)) - vect2*rve_stretch(g(1)) ;
        bnd_nodesf(:,g(2)) = bnd_nodes(:,g(2)) + vect2*rve_stretch(g(2)) ;
    elseif p == 5
        nodesf(:,g(1)) = nodes(:,g(1)) - vect1*rve_stretch(g(1)) ;
        bnd_nodesf(:,g(1)) = bnd_nodes(:,g(1)) - vect2*rve_stretch(g(1)) ;
    elseif p == 6
        nodesf(:,g(1)) = nodes(:,g(1)) - vect1*rve_stretch(g(1)) ;
        nodesf(:,g(2)) = nodes(:,g(2)) - vect1*rve_stretch(g(2)) ;
        bnd_nodesf(:,g(1)) = bnd_nodes(:,g(1)) - vect2*rve_stretch(g(1)) ;
        bnd_nodesf(:,g(2)) = bnd_nodes(:,g(2)) - vect2*rve_stretch(g(2)) ;
    elseif p == 7
        nodesf(:,g(2)) =  nodes(:,g(2)) - vect1*rve_stretch(g(2)) ;
        bnd_nodesf(:,g(2)) = bnd_nodes(:,g(2)) - vect2*rve_stretch(g(2)) ;
    elseif p == 8
        nodesf(:,g(1)) = nodes(:,g(1)) + vect1*rve_stretch(g(1)) ;
        nodesf(:,g(2)) = nodes(:,g(2)) - vect1*rve_stretch(g(2)) ;
        bnd_nodesf(:,g(1)) = bnd_nodes(:,g(1)) + vect2*rve_stretch(g(1)) ;
        bnd_nodesf(:,g(2)) = bnd_nodes(:,g(2)) - vect2*rve_stretch(g(2)) ;
    end
    
    plot3(nodesf(:,1) , nodesf(:,2), nodesf(:,3),'o', 'LineWidth',0.2, 'MarkerEdgeColor','k','MarkerFaceColor',[0 0.5 0.8],'MarkerSize',5);
    
    m = 1 ;
    
    for n = 1 : length(fibers) % count rows
        
        node1 = nodesf(fibers(n,1),:) ;
        
        for i = 1:3
            for k = 1:abs(fibers(n,2+i))
                node2 = bnd_nodesf(m,:) ;
                
                x1(1) = node1(1); % node 1 x coord
                y1(1) = node1(2); % node 1 y coord
                z1(1) = node1(3); % node 1 z coord
                
                x1(2) = node2(1); % node 2 x coord
                y1(2) = node2(2); % node 2 y coord
                z1(2) = node2(3); % node 2 z coord
                
                if fib_type(n) == 1
                    C = [0.8 0.8 0.4];
                    lw = 1.5;
                elseif fib_type(n) == 2
                    C = [0 0.2 0];
                    lw=2;
                elseif fib_type(n) == 3
                    C = [0.7 0 0];
                    lw=1.5;
                end
                
                plot3(x1, y1, z1, 'Color', C,'LineWidth',lw);
                
                m = m+1 ;
                
                node1 = bnd_nodesf(m,:) ;
                
                m = m+1 ;
            end
        end
        
        node2 = nodesf(fibers(n,2),:) ;
        
        x1(1) = node1(1); % node 1 x coord
        y1(1) = node1(2); % node 1 y coord
        z1(1) = node1(3); % node 1 z coord
        
        x1(2) = node2(1); % node 2 x coord
        y1(2) = node2(2); % node 2 y coord
        z1(2) = node2(3); % node 2 z coord
        
        if fib_type(n) == 1
            C = [0.8 0.8 0.4];
            lw = 1.5;
        elseif fib_type(n) == 2
            C = [0 0.2 0];
            lw = 2;
        elseif fib_type(n) == 3
            C = [0.7 0 0];
            lw = 1.5;
        end
        
        plot3(x1, y1, z1, 'Color', C, 'LineWidth',lw);
    end
end

%% plotting outline around the central RVE
row1x = [-0.5,(-0.5+rve_stretch(1))];
row2x = (-0.5+rve_stretch(1)).*ones(length(row1x),1) ;
lowrow = [-0.5,-0.5] ;
row1y = [-0.5,(-0.5+rve_stretch(2))];
row2y = (-0.5+rve_stretch(2)).*ones(length(row1y),1) ;
row1z = [-0.5,(-0.5+rve_stretch(3))];
row2z = (-0.5+rve_stretch(3)).*ones(length(row1z),1) ;

lw = 2;

plot3(row1x, lowrow, lowrow,'k--','LineWidth', lw)
plot3(row1x, row2y, lowrow,'k--','LineWidth', lw)
plot3(row1x, lowrow, row2z,'k--','LineWidth', lw)
plot3(row1x, row2y, row2z,'k--','LineWidth', lw)
plot3(lowrow, row1y, lowrow,'k--','LineWidth', lw)
plot3(row2x, row1y, lowrow,'k--','LineWidth', lw)
plot3(lowrow, row1y, row2z,'k--','LineWidth', lw)
plot3(row2x, row1y, row2z,'k--','LineWidth', lw)
plot3(lowrow, lowrow, row1z,'k--','LineWidth', lw)
plot3(lowrow, row2y, row1z,'k--','LineWidth', lw)
plot3(row2x, lowrow, row1z,'k--','LineWidth', lw)
plot3(row2x, row2y, row1z,'k--','LineWidth', lw)

%% formatting and labelling plot
set(gcf, 'color', 'white');
axis square;
xlabel('Circumferential')'; ylabel('Axial'); zlabel('Radial');

set(gca, 'FontSize',16)

end
