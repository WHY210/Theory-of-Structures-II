%%% 結構二加分題-勁度法-例17.2 %%% 
%%% 土木13B-吳巽言-109612054  %%%


%% 參數輸入區 (假設E=A=I=L=1,其他同題目敘述)

% member 1 
E1 = 1;
A1 = 1;
I1 = 1;
L1 = 1;
x1 = acos(0.6);
y1 = acos(0.8);

% member 2 
E2 = 1;
A2 = 1;
I2 = 1;
L2 = 1;
x2 = acos(1);
y2 = acos(0);

% member 3           
E3 = 1;
A3 = 1;
I3 = 1;
L3 = 1;
x3 = acos(0);
y3 = acos(1);

% Force
  % [Q1, Q2, Q3, Q4]
Q = [ 0;  0;  0;  1];

% code
    % [Nx, Ny, Nz, Fx, Fy, Fz]
mg1 = [ 0,  0,  0,  4,  4,  1];
mg2 = [ 4,  4,  2,  4,  0,  3];
mg3 = [ 0,  0,  0,  4,  0,  3];
    % [Nx', Ny', Nz', Fx', Fy', Fz']
ml1 = [  0,   0,   0,   0,   4,   1];
ml2 = [  4,   4,   2,   4,   0,   3];
ml3 = [  0,   0,   0,   0,   4,   3];

% lateral displacement
    % [Nx,  Ny, Nz, Fx,  Fy, Fz]
lg1 = [ 1,   1,  1,  1,-3/4,  1];
lg2 = [ 1,-3/4,  1,  1,   1,  1];
lg3 = [ 1,   1,  1,  1,   1,  1];
    % [Nx', Ny', Nz', Fx',  Fy', Fz']
ll1 = [ 1,    1,   1,   1, -5/4,   1];
ll2 = [ 1, -3/4,   1,   1,    1,   1];
ll3 = [ 1,    1,   1,   1,   -1,   1];



%% compute

% member stiffness matrices 
[kg1,kl1] = StiffnessMatric(E1, A1, I1, L1, x1, y1);
[kg2,kl2] = StiffnessMatric(E2, A2, I2, L2, x2, y2);
[kg3,kl3] = StiffnessMatric(E3, A3, I3, L3, x3, y3);

% structure stiffness matrices
K = GlobalStiffnessMatrix(kg1,kg2,kg3,mg1,mg2,mg3,lg1,lg2,lg3);

% Q=KD => Displacement:[D1,D2,D3,D4]
D = inv(K)*Q;

% q=k'd => force in member:[qNx', qNy', qNz', qFx',qFy', qFz']
q1 = ForceInMember(ml1,ll1,D,kl1);
q2 = ForceInMember(ml2,ll2,D,kl2);
q3 = ForceInMember(ml3,ll3,D,kl3);



%% results
disp("K = "); disp(K);
disp("D = "); disp(D);
disp("q1 = "); disp(q1);
disp("q2 = "); disp(q2);
disp("q3 = "); disp(q3);



%% functions

function [kg,kl] = StiffnessMatric(E, A, I, L, x, y)
    % member stiffness matrix in local coordinates 
    kl = [ A*E/L,           0,          0, -A*E/L,           0,          0; 
               0,  12*E*I/L^3,  6*E*I/L^2,      0, -12*E*I/L^3,  6*E*I/L^2;
               0,   6*E*I/L^2,   4*E*I/L,       0,  -6*E*I/L^2,    2*E*I/L;
          -A*E/L,           0,          0,  A*E/L,           0,          0;
               0, -12*E*I/L^3, -6*E*I/L^2,      0,    12*E*I/L, -6*E*I/L^2;
               0,   6*E*I/L^2,    2*E*I/L,      0,  -6*E*I/L^2,    4*E*I/L ];

    % displacement transformation matrix
    T = [  cos(x), cos(y), 0,       0,      0, 0;
          -cos(y), cos(x), 0,       0,      0, 0;
                0,      0, 1,       0,      0, 0;
                0,      0, 0,  cos(x), cos(y), 0;
                0,      0, 0, -cos(y), cos(x), 0;
                0,      0, 0,       0,      0, 1 ];

    % force transformation matrix
    transT = transpose(T);
    
    % member stiffness matrix in global coordinates
    kg = transT * kl * T;
end

% structure stiffness metrix K
function K = GlobalStiffnessMatrix(kg1,kg2,kg3,m1,m2,m3,lg1,lg2,lg3)
    i=1;
    j=1;
    K = zeros(4,4);
    while i<=4
        while j<=4         
            i1 = find(m1 == i); j1 = find(m1 == j); 
            i2 = find(m2 == i); j2 = find(m2 == j); 
            i3 = find(m3 == i); j3 = find(m3 == j); 
           
            kij1 = lg1(i1)*kg1(i1,j1)*transpose(lg1(j1)); 
            kij2 = lg2(i2)*kg2(i2,j2)*transpose(lg2(j2)); 
            kij3 = lg3(i3)*kg3(i3,j3)*transpose(lg3(j3));
            K(i,j) = sum(transpose(sum(kij1))) ...
                   + sum(transpose(sum(kij2))) ...
                   + sum(transpose(sum(kij3)));
            j = j+1;     
        end
        i = i+1;
        j = 1;
    end
end
        
function q = ForceInMember(ml,ll,D,kl)
  
    % member local displacement
    d = [];
    for a = 1:6
        ad = ml(a); 
        if ad == 0
            d(a) = 0;
        else
            d(a) = D(ad)*ll(a);
        end     
    end

    % q=k'd => member force q
    q = kl*transpose(d); 
end