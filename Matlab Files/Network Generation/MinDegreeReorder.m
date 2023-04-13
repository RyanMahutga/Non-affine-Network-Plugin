function [ nodes_n, fibers_n, stab ] = MinDegreeReorder( nodes, fibers, stab_nodes )
%MinDegreeReorder.m reorders the nodes so that the Jacobian has minimal bandwidth

%   INPUTS: nodes - Nx3 matrix of nodal coordinates'
%            fibers - Mx5 matrix of nodal connectivity and crossing conditions
             
%   OUTPUTS: nodes_n - Nx3 matrix of reordered nodes
%            fibers_n - Mx5 matrix of fibers including new node numbers in
%            (m,1) and (m,2)

% NOTE: This algorithm doesn't necessarily work properly for fiber networks
% with multiple fibers connecting the same two nodes. (i.e. the nodal
% connectivity changes for this case, and that connectivity is only counted
% as 1 in this algorithm)

% Created by: Ryan Mahutga - Barocas Research Group - University of Minnesota
% Date Created: 02-12-18

N=length(nodes) ;
M=length(fibers) ;

stab = stab_nodes;

nodes_n = zeros(N,3) ;

B = zeros(N,N) ;

for i=1:N
    B(i,i) = 1 ;
end

for j = 1:M
   idx1 = fibers(j,1) ;
   idx2 = fibers(j,2) ;
   
   B(idx1,idx2) = 1 ;
   B(idx2,idx1) = 1 ;
end

p = symamd(B) ;

fibers_n = fibers;

for k = 1:N
    nodes_n(k,:) = nodes(p(k),:) ;
    
    if stab_nodes(1) == p(k)
        stab(1) = k;
    elseif stab_nodes(2) == p(k)
        stab(2) = k;
    elseif stab_nodes(3) == p(k)
        stab(3) = k;
    end
    
    fibers_n(fibers(:,1)==p(k),1)=k ;
    fibers_n(fibers(:,2)==p(k),2)=k ;
end

% subplot(2,2,1), spy(B), title('Before')
% subplot(2,2,2), spy(B(p,p)), title('After')

end

