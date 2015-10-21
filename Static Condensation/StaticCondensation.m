%Static Condensation Method

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
    
    %Transfomation Matrix (Tr) (square matrix will be used)
    Tr = zeros(numberNodes);
    Tr(1,1) = 1;
    Tr(2,2) = 1;
    Tr(3,3) = 1;
    Tr(4,3) = 1;
    Tr(4,5) = 1;
    Tr(5,5) = 1;
    Tr(6,6) = 1;
    Tr(7,7) = 1;
    UF = U;
    KF = transpose(Tr)*K*Tr;
    KF(4,4) = 1;
    FF = transpose(Tr)*F;
    TF = transpose(Tr)*T*Tr;
    TF(4,4) = 1;
    
    
    %Begin Newton Method
    KFred = KF(2:(numberNodes), 2:(numberNodes));
    TFred = TF(2:(numberNodes), 2:(numberNodes));
    FFred = FF(2:(numberNodes));
    UFred = UF(2:(numberNodes));
    
    %Residual
    r = KFred*UFred - FFred; 
    
    U0 = 0;
    UF_new = UFred - TFred\r;
    condition = norm(UF_new - UFred)/norm(UF_new)
    U = [U0;UF_new];
    U(4,1) = U(3,1) + U(5,1);
    iter(count) = count;
    normres(count) = norm(r);
    normu(count) = condition;
end


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
