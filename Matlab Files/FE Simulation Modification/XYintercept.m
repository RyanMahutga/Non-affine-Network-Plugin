function [ x,y ] = XYintercept( node1, node2)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

lend = norm(node2-node1); % length of fiber
d = (node2-node1)./lend; % direction vector

x = d(1)/d(3) * (0-node1(3)) + node1(1);
y = d(2)/d(3) * (0-node1(3)) + node1(2);

end

