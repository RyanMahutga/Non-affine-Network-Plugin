function [ fib_rads, init_lens ] = networkFeatures( nodes, fibers, mu, sigma)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

prad = makedist('Normal','mu', mu, 'sigma', sigma) ;

for d = 1:length(fibers)
    node1 = nodes(fibers(d,1),:) ;
    node2 = nodes(fibers(d,2),:) + fibers(d,3:5);
    
    init_lens(d) = norm(node2-node1) ;
    
    fib_rads(d)=mu;
%     while fib_rads(d) <=0
%         fib_rads(d) = random(prad,1,1) ;
%     end
end

end

