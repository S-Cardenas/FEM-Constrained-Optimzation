%Penalty Method

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

%for structure
%displacements: displacement vector (U)
U=zeros(numberNodes,1);

condition = 1;
condition_res = 1;
count = 0;

%Creation of Element Stiffness Matrix (K) and Tangent Stiffness Matrix (T)
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
    
    %Assembly of Global Stiffness Matrices
    for e=1:numberElements
        node1 = elementNodes(e,1);
        node2 = elementNodes(e,2);
        Ke = elementStiffness(U, node1, node2);
        Te = elementTangentStiffness(U, node1, node2);
        K(node1:node2, node1:node2) = K(node1:node2, node1:node2) + Ke;
        T(node1:node2, node1:node2) = T(node1:node2, node1:node2) + Te;
    end
    
    %Creation of Tangent Constraint Matrix (BT)
    B = zeros(2,7);
    B(1,1) = 1;
    B(2,3) = 1;
    B(2,4) = -1;
    B(2,5) = 1;
    
    %Assign Penalty Parameter (Beta)
    Beta = 1e11;
    
    %Assign Vector of Constants (Q)
    Q = zeros(2,1);
    
    %Compile Individual Stiffness Matrices, Displacement vector, and Force 
    %Vector into Final Stiffness Matrix(KF), Final Displacement Vector (UF)
    %, and Final Force Vector (FF).
    KF = K + Beta*transpose(B)*B;
    UF = U;
    FF = F + Beta*transpose(B)*Q;
    
    %Compile Tangent Final Stiffness Matrix (TF)
    TF = T + Beta*transpose(B)*B;
    
    %Begin Newton Method
    KFred = KF(1:numberNodes, 1:numberNodes);
    TFred = TF(1:numberNodes, 1:numberNodes);
    FFred = FF(1:numberNodes);
    UFred = UF(1:numberNodes);
    
    r = KFred*UFred - FFred; %Residual
    
    UF_new = UFred - TFred\r;
    condition = norm(UF_new - UFred)/norm(UF_new);
    U = [UF_new];
    iter(count) = count;
    normres(count) = norm(r);
    normu(count) = condition;
end


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
