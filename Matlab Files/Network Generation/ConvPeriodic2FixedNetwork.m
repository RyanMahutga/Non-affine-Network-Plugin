function ConvPeriodic2FixedNetwork( nodes, bnd_nodes, fibers, fib_type )
%ConvPeriodic2FixedNetwork converts a period network to a fixed network for
% use comparing the periodic function to the standard fixed cauchy solution
%   INPUTS: nodes-Nx3 matrix of nodal coordinates (x,y,z)
%           fibers-Mx5 matrix of node numbers on the ends of fibers and
%           fiber crossings
%   OUTPUTS: NONE

% Ryan Mahutga - Barocas Research Group - University of Minnesota
% Date Created: 1-23-18

fibers_n1 = [];
c = length(nodes) ;

m = 1 ;
fibtype=[];
for n = 1 : length(fibers) % count rows
        
        node1 = nodes(fibers(n,1),:) ;        
        
        fibers_n1 = [fibers_n1, fibers(n,1)] ;
        
        for i = 1:3
            for k = 1:abs(fibers(n,2+i))                
                node2 = bnd_nodes(m,:) ;                                
                nodes = [nodes; node2] ;
                fibers_n1 = [fibers_n1,c+m];   
                m = m+1 ;
                
                node1 = bnd_nodes(m,:) ;
                nodes = [nodes; node1] ;
                fibers_n1 = [fibers_n1, c+m];
                m = m+1 ;
                fibtype = [fibtype, fib_type(n)];
            end
        end
        
        fibtype = [fibtype, fib_type(n)];
        fibers_n1 = [fibers_n1, fibers(n,2)] ;       
       
end

fibers_n = zeros(length(fibers_n1)/2,2) ;
fibers_n(:,1) = fibers_n1(1:2:end) ;
fibers_n(:,2) = fibers_n1(2:2:end) ;

init_lens = zeros(1,length(fibers_n)) ;

for n = 1:length(fibers_n)
   init_lens(1,n) = pdist([nodes(fibers_n(n,1),:) ; nodes(fibers_n(n,2),:)]) ; 
end

fib_rads = ones(length(fibers_n),1) ;
fib_rads = 30e-9.*fib_rads;

fibers = [];
fibers = fibers_n ;

phi=0.15;

put_net(nodes, fibers, fibtype', init_lens', fib_rads, phi)

mkdir('FixedNetwork')
old=cd('FixedNetwork');

save('nodes.mat','nodes') 
save('fib_rads.mat','fib_rads') 
save('fibers.mat','fibers')
save('fibtype.mat', 'fibtype')
save('init_lens.mat', 'init_lens')

cd(old)

end

