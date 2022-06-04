 clc 
clearvars
% Lo anterior para empezar cualquier codigo limpiando las variables y la
% ventana de comandos
% Cedula jmotoar: 1008547364
% cedula abuitragol: 1053869374
% Auxiliares para entender mejor el codigo
X = 1; Y = 2; 
% P = ceil((4+4)/2);

% Numero de nodos, elementos y grados de libertad
nnod = 3;    nelem = 2;   ngdl = nnod*2;

%% ------- (1) DATOS DE ENTRADA --------
 
% ***** DATOS DE LOS NODOS ******

% Coordenadas de los nodos: fila=nodo_i, coll=Xi(m), col2=Uy

coord = [0  0
         7  0
         12 0];

% Grados de libertad fila=nodo_i, coll=Uxi, col2=Uy

% gdl = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12];

%ciclo for 
gdl = zeros(nnod,2);  %se separa la memoria, dos columnas por Ux Uy
for i=1:nnod
    gdl(i,:) = [(2*i-1) 2*i];
end

%% ******** DATOS DE LOS ELEMENTOS *******

b = 0.1; h = 0.2;
A = b*h;      % Area de los elementos en m^2
E = 200000000;  % Modulo de elasticidad en Kn/m
I = (b*(h^3))/12;
% Elementos: fila = elemento_e, colum1 = nodo_i, colum2 = nodo_j 

elem = [1 2
        2 3];

% longitudes y angulo de orientacion de cada elemento
% Se separa la memoria:
L =  zeros(nelem,1);

for e = 1:nelem

    % Se extraen los nodos i y j de cada elemento
    nodo_i = elem(e,1);
    nodo_j = elem(e,2);

    % Se obtienen las coordenadas de los nodos i y j
    Xi = coord(nodo_i,X);   
    Yi = coord(nodo_i,Y);   
    Xj = coord(nodo_j,X);
    Yj = coord(nodo_j,Y);

    % Se calculan las longitudes en la eq(29)
    L(e) = sqrt((Xj-Xi)^2 + (Yj-Yi)^2);

end

% MATRIZ LaG: cada fila es un elemento e
LaG = zeros(nelem,4);
for e = 1:nelem
    LaG(e,:) = LaG(e,:) + [(elem(e,1)*2)-1, elem(e,1)*2, (elem(e,2)*2)-1, elem(e,2)*2 ];
end

% ******** DATOS DE LOS APOYOS Y CARGAS ********

% Grados de libertad restringidos(conocidos)
a = [1 2 3 5 6];

% Despazamientos gdl restringidos
Da = [0 0 0 0 0 ]';

% Fuerzas en los gdl restringidos (apoyos) 
Faa = [0 0 0 0 0]';

% Grados de libertad no restringidos
b = 4;

% Fuerzas en los gdl no restringidos (fuerzas nodales externas)
Pb = 0;

% Vectores de fuerzas de empotramiento en los elementos

f1 = [70 81.67 70 -81.67]';
f2 = [14.08 19.20 25.92 -28.80]';

fe = [f1 f2];
%% ------ (2) PROCESO DE CALCULO --------

% Se separa la memoria
ke = cell(nelem,1);
K = zeros(ngdl, ngdl);
F = zeros(ngdl , 1);

% ****** ENSAMBLAJE MATRICIAL ******

% Con ciclo for para optimizar el proceso ya que es un proceso repetitivo
% Matriz de rigidez ke en sistema local eq(17)
for e = 1:nelem
ke{e} = [12*E*I/L(e)^3 6*E*I/L(e)^2 -12*E*I/L(e)^3 6*E*I/L(e)^2
        -6*E*I/L(e)^2 4*E*I/L(e) -6*E*I/L(e)^2 2*E*I/L(e)
        -12*E*I/L(e)^3 -6*E*I/L(e)^2 12*E*I/L(e)^3 -6*E*I/L(e)^2
         6*E*I/L(e)^2 2*E*I/L(e) -6*E*I/L(e)^2 4*E*I/L(e)];

% Matriz de rigidez Ke en sistema global eq(25)
Ke = ke{e};
% Grados de libertad globales del elemento e
gdl_e = LaG(e,:);
% Se suma las contribuciones ke a la matriz K
K(gdl_e, gdl_e) = K(gdl_e, gdl_e) + Ke;
% Ensamblaje vcetor de fuerzas de empotramiento 
F(gdl_e,1) = F(gdl_e,1) + fe(:,e);
end



% Matriz de rigidez total K [KN/m]
%disp(k);
%sky(k);

% ***** PARTICION DE SUBMATRICES *****
%
% ( Pa )     ( Kaa | Kab )  (Da)   (Fa)
%-------  =  ------|------  ---- + ----
% ( Pb )     ( Kba | Kbb )  (Db)   (Fb)

Kaa = K(a,a);
Kab = K(a,b);
Kba = K(b,a);  
Kbb = K(b,b);
Fa = F(a) - Faa;
Fb = F(b);

% Se calculan los desplazamientos desconocidos eq(7) y las reacciones
% eq(18)
Db = Kbb\(Pb - Kba*Da -Fb);
Pa = Kaa*Da  + Kab*Db + Fa;
% se arman los vectores D y P
D = zeros(ngdl,1);
D(a) = Da;
D(b) = Db;
P = zeros(ngdl,1); 
P(a) = Pa; 
P(b) = Pb; 

% Fuerzas internas de cada elemento Eq(36)
pe = zeros(4,nelem);
V = zeros(nelem,2); % Cortantes
M = zeros(nelem,2); % Momentos
for e = 1:nelem
    De = D(LaG(e,:));

    pe(:,e) = ke{e}*De + fe(:,e);
    %          Vi        Vj
    V(e,:) = [pe(1,e) pe(3,e)];
    %          Mi        Mj
    M(e,:) = [pe(2,e) pe(4,e)];
end

%% ************************ MUESTRA DE RESULTADOS ************************

format shortG
disp('-------- RESULTADOS --------')
disp('')

% MOSTRAR CORTANTES
disp('Estas son las Fuerzas Cortantes en los apoyos')
disp(V)
disp('')

%MOSTRAR MOMENTOS
disp('Estos son los Los momentos flectores de cada nodo')
disp(M)
disp('')

% %MOSTRAR FUERZAS AXIALES DE CADA ELEMENTO
% disp('Estas son las fuerzas axiales de cada elemento')
% disp(N)