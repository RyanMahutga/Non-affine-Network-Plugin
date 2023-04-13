function [ x,y ] = XYintercept( node1, node2, dir)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

lend = norm(node2-node1); % length of fiber
d = (node2-node1)./lend; % direction vector

if dir == 3
    x = d(1)/d(3) * (0-node1(3)) + node1(1);
    y = d(2)/d(3) * (0-node1(3)) + node1(2);
elseif dir==2
    x = d(2)/d(1) * (0-node1(1)) + node1(2);
    y = d(3)/d(1) * (0-node1(1)) + node1(3);
elseif dir==1
    x = d(1)/d(2) * (0-node1(2)) + node1(1);
    y = d(3)/d(2) * (0-node1(2)) + node1(3);
end

end

