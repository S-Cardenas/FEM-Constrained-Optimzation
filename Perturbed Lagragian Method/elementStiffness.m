function [Ke] = elementStiffness(U, node1, node2)
u1 = U(node1);
u2 = U(node2);
Ke(1,1) = 2+u2-u1;
Ke(1,2) = -2-u2+u1;
Ke(2,1) = -2-u2+u1;
Ke(2,2) = 2+u2-u1;
