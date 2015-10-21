%Lagrange Multiplier Method

%clear memory
clear all

%element length
L=1; 

%Tolerance
tol=1.0e-3;

%numberElements: number of Elements
numberElements=6;

%numberNodes: number of Nodes
numberNodes=7;

%elementNodes: connections at elements
ii=1:numberElements;
elementNodes(:,1)=ii;
elementNodes(:,2)=ii+1;

%displacements: displacement vector (U)
U=zeros(numberNodes,1);

%Number of Legrange Multipliers
NumLegrange = 2;
    
%Creation of Largrange Multiplier Matrix (L)
L = zeros(NumLegrange,1);

condition = 1;
condition_res = 1;
count = 0;

while (condition > tol)
    count = count+1
    
    %for structure
        %force: force vector (F)
        %K: Global stiffness matrix
        %T: Global tangent stiffness matrix
    F=zeros(numberNodes,1);
    F(2) = 1;
    F(7) = -1;
    K=zeros(numberNodes);
    T=zeros(numberNodes);
    
    %Assembly of Global Stiffness Matrix (K) and Globl Tangent Stiffness
    %Matrix (T)
    for e=1:numberElements
        node1 = elementNodes(e,1);
        node2 = elementNodes(e,2);
        Ke = elementStiffness(U, node1, node2);
        Te = elementTangentStiffness(U, node1, node2);
        K(node1:node2, node1:node2) = K(node1:node2, node1:node2) + Ke;
        T(node1:node2, node1:node2) = T(node1:node2, node1:node2) + Te;
    end
    
    %Creation of Constraint Matrix (B)
    B = zeros(NumLegrange,numberNodes);
    B(1,1) = 1;
    B(2,3) = 1;
    B(2,4) = -1;
    B(2,5) = 1;
    
    %Creation of Tangent Constraint Matrix (BT)
    BT = zeros(NumLegrange,numberNodes);
    BT(1,1) = 1;
    BT(2,3) = 1;
    BT(2,4) = -1;
    BT(2,5) = 1;
    
    %Creation of zero marix (Z)
    Z = zeros(NumLegrange,NumLegrange);
    
    %Creation of Tangent zero marix (ZT)
    ZT = zeros(NumLegrange,NumLegrange);
    
    %Creation of Largrange Equivalent Force Matrix (N)
    N = zeros(NumLegrange,1);
    
    %Compile Individual Stiffness Matrices, Displacement vector, and Force 
    %Vector into Final Stiffness Matrix(KF), Final Displacement Vector (UF)
    %, and Final Force Vector (FF).
    KF = [K, transpose(B); B, Z];
    UF = [U; L];
    FF = [F; N];
    
    %Compile Tangent Final Stiffness Matrix (TF)
    TF = [T, transpose(BT); BT, ZT];
    
    %Begin Newton Method
    KFred = KF(1:(numberNodes+NumLegrange), 1:(numberNodes+NumLegrange));
    TFred = TF(1:(numberNodes+NumLegrange), 1:(numberNodes+NumLegrange));
    FFred = FF(1:(numberNodes+NumLegrange));
    UFred = UF(1:(numberNodes+NumLegrange));
    
    %Residual
    r = KFred*UFred - FFred; 
    
    UF_new = UFred - TFred\r;
    condition = norm(UF_new - UFred)/norm(UF_new)
    U = UF_new(1:numberNodes);
    L = UF_new(numberNodes+1:numberNodes+NumLegrange);
    iter(count) = count;
    normres(count) = norm(r);
    normu(count) = condition;
end


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
