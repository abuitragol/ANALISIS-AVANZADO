E = 5;
A = 3;
L = 2;
I = 4;
types = [1 1 1 1];
ket = cell {4,1};
tipos_carga = {puntual momento_par distri distrit axial axialdistri};
% tipo 1
ket{1} = [E*A/L,           0,            0, -E*A/L,           0,             0
            0,  (12*E*I)/L^3,  (6*E*I)/L^2,     0, -(12*E*I)/L^3,   (6*E*I)/L^2
            0,   (6*E*I)/L^2,    (4*E*I)/L,     0,  -(6*E*I)/L^2,     (2*E*I)/L
        -E*A/L,            0,            0, E*A/L,            0,             0
            0, -(12*E*I)/L^3, -(6*E*I)/L^2,     0,  (12*E*I)/L^3,  -(6*E*I)/L^2
            0,   (6*E*I)/L^2,    (2*E*I)/L,     0,  -(6*E*I)/L^2,    (4*E*I)/L];
% tipo 2
ket{2} = [E*A/L,            0,           0, -E*A/L,            0,            0
             0,   (3*E*I)/L^3,           0,      0,  -(3*E*I)/L^3,   (3*E*I)/L^2
             0,             0,           0,      0,             0,           0
         -E*A/L,            0,           0,  E*A/L,            0,            0
             0,  -(3*E*I)/L^3,           0,      0,   (3*E*I)/L^3,  -(3*E*I)/L^2
             0,   (3*E*I)/L^2,           0,      0,  -(3*E*I)/L^2,    (3*E*I)/L];
% tipo 3
ket{3} = [E*A/L,            0,            0, -E*A/L,            0,           0
             0,   (3*E*I)/L^3,  (3*E*I)/L^2,      0, -(3*E*I)/L^3,           0 
             0,   (3*E*I)/L^2,    (3*E*I)/L,      0, -(3*E*I)/L^2,           0
         -E*A/L,            0,            0,  E*A/L,            0,           0
             0,  -(3*E*I)/L^3, -(3*E*I)/L^2,      0,  (3*E*I)/L^3,           0
             0,             0,            0,      0,            0,           0];
% tipo 4
ket{4} = [   E*A/L,            0,           0, -E*A/L,            0,           0
                 0,            0,           0,      0,            0,           0
                 0,            0,           0,      0,            0,           0
            -E*A/L,            0,           0,  E*A/L,            0,           0
                 0,            0,           0,      0,            0,           0
                 0,            0,           0,      0,            0,           0];

%vectores de fuerzasde empotramiento
fet = cell{4,1};
fet{1} = [f1
          f2
          f3
          f4
          f5
          f6];

fet{2} = [f1
          f2-((3/(2*L(e)))*f3)
          0
          f4
          f5+((3/(2*L(e)))*f3)
          f6-((1/2)*f3)];
fet{3} = [f1
          f2-((3/(2*L(e)))*f6)
          f3-((1/2)*f6)
          f4
          f5+((3/(2*L(e)))*f6)
          0];
fet{4} = [f1
          f2-((1/L(e)*(f3+f6)))
          0
          f4
          f5-((1/L(e)*(f3+f6)))
          0];
%matriz de transformacion de porticos

T = [   eta(e),     mu(e),           0,      0,            0,           0
        -mu(e),    eta(e),           0,      0,            0,           0
          0,            0,           1,      0,            0,           0
          0,            0,           0,    eta(e),        mu(e),        0
          0,            0,           0,    -mu(e),       eta(e),        0
          0,            0,           0,      0,            0,           1];

%creacion de memoria de vector de fuerzas 
Fet = cell(nelem,1);
for e = 1:nelem
    Fet{e} = fet{types,1};
end