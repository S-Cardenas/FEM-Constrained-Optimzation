function [Te] = elementTangentStiffness(U, node1, node2)
u1 = U(node1);
u2 = U(node2);
Te(1,1) = 2+2*u2-2*u1;
Te(1,2) = -2-2*u2+2*u1;
Te(2,1) = -2-2*u2+2*u1;
Te(2,2) = 2+2*u2-2*u1;