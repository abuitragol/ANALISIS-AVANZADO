clc 
clearvars
% Lo anterior para empezar cualquier codigo limpiando las variables y la
% ventana de comandos
% Cedula jmotoar: 1008547364
% cedula abuitragol: 1053869374
% Auxiliares para entender mejor el codigo
X = 1; Y = 2; P = ceil((4+4)/2);

% Numero de nodos, elementos y grados de libertad
nnod = 15;    nelem = 27;   ngdl = 30;

%% ------- (1) DATOS DE ENTRADA --------
 
% ***** DATOS DE LOS NODOS ******

% Coordenadas de los nodos: fila=nodo_i, coll=Xi(m), col2=Uy

coord = readmatrix("COORDENADAS CERCHAS.xlsx","Range",'B2:C16');

% Grados de libertad fila=nodo_i, coll=Uxi, col2=Uy

% gdl = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12];

%ciclo for 
gdl = zeros(nnod,2);  %se separa la memoria, dos columnas por Ux Uy
for i=1:nnod
    gdl(i,:) = [(2*i-1) 2*i];
end

%% ******** DATOS DE LOS ELEMENTOS *******

A = 0.001495;      % Area de los elementos en m^2
E = 200000000;  % Modulo de elasticidad en Kn/m

% Elementos: fila = elemento_e, colum1 = nodo_i, colum2 = nodo_j 

elem = readmatrix("COORDENADAS CERCHAS.xlsx","Range",'G2:H28');

% longitudes y angulo de orientacion de cada elemento
% Se separa la memoria:
L =  zeros(nelem,1);
eta = zeros(nelem,1);
mu = zeros(nelem,1);

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

    % Y los cosenos y senos directores con las eq(27) y eq(28)
    eta(e) = (Xj-Xi)/L(e);
    mu(e) = (Yj-Yi)/L(e);
end

% MATRIZ LaG: cada fila es un elemento e

LaG = readmatrix("COORDENADAS CERCHAS.xlsx","Range",'O2:R28');

% ******** DATOS DE LOS APOYOS Y CARGAS ********

% Grados de libertad restringidos(conocidos)
a = [1 2 11 12];

% Despazamientos gdl restringidos
Da = [0 0 0 0]';

% Grados de libertad no restringidos
b = [3 4 5 6 7 8 9 10 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 ];

% Fuerzas en los gdl no restringidos (fuerzas nodales externas)
Pb = [0 0 0 0 0 0 0 0 0 -20*P -(20*P)*cosd(69.31) -(20*P)*sind(69.31) 0 -30*P 0 0 0 0 0 -30*P (40*P)*cosd(69.31) -(40*P)*sind(69.31) 0 -40*P 0 -50*P ]';

% Fuerzas en los gdl de los apoyos 
Fa = [0 -10*P 0 -10*P]';

%% ------ (2) PROCESO DE CALCULO --------

% Se separa la memoria
K = zeros(ngdl, ngdl);

% ****** ENSAMBLAJE MATRICIAL ******

% Con ciclo for para optimizar el proceso ya que es un proceso repetitivo
% Matriz de rigidez ke en sistema local eq(17)
for e = 1:nelem
ke = [A*E/L(e) 0 -A*E/L(e) 0;0 0 0 0; -A*E/L(e) 0 A*E/L(e) 0; 0 0 0 0];

% MatriZ de transformacion de cada elemento eq(22)

T = [eta(e) mu(e) 0 0; -mu(e) eta(e) 0 0; 0 0 eta(e) mu(e); 0 0 -mu(e) eta(e)];

% Matriz de rigidez Ke en sistema global eq(25)
Ke = T'*ke*T;
% Grados de libertad globales del elemento e
gdl_e = LaG(e,:);
% Se suma las contribuciones ke a la matriz K
K(gdl_e, gdl_e) = K(gdl_e, gdl_e) + Ke;
end



% Matriz de rigidez total K [KN/m]

% ***** participacion de submatrices *****

Kaa = K(a,a);
Kab = K(a,b);
Kba = K(b,a);  
Kbb = K(b,b);

% Se calculan los desplazamientos desconocidos eq(7) y las reacciones
% eq(18)
Db = Kbb\(Pb - Kba*Da);
Pa = (Kaa*Da  + Kab*Db)-Fa;
% se arman los vectores D y P
D = zeros(ngdl,1);
D(a) = Da;  D(b) = Db;
P = zeros(ngdl,1); 
P(a) = Pa; 
P(b) = Pb;

% ********** CALCULO DE FUERZAS AXIALES **********

% Fuerzas axiales internas de cada elemento eq (30)
N = zeros(nelem,1);
for e = 1:nelem
    %  Extraemos los elementos del vector de desplazamientos
    De = D(LaG(e,:));
% Se calcula respuesta interna de cada elemento              
     N(e) = ((E*A)/L(e))*[-eta(e) -mu(e) eta(e) mu(e)]*De;
% Se agrega la respuesta de cada elemento al vector de fuerzas internas 
     N(e,1)= N(e);
end
%% ****************** GRAFICAS ******************
figure %comando para empezar una figura
for e = 1:nelem
    hold on % hold on es el comando para que mantenga la grafica del punto anterior y sobre este se sobreescriba el siguiente
    m = elem(e, X); % Nodo i
    n = elem(e, Y); % Nodo j

    x_i = coord(m, X); % Extraigo la coordenada en X del nodo i
    y_i = coord(m, Y); % Extraigo la coordenada en Y del nodo i 
    x_j = coord(n, X); % Extraigo la coordenada en X del nodo j
    y_j = coord(n, Y); % Extraigo la coordenada en Y del nodo j
    plot([x_i, x_j], [y_i, y_j], 'k') % Grafico las coordenadas de cada nodo y se unen
end

D = D*3; % Magnificacion X3 de los desplazamientos para efectos visuales en la animacion

for paso = 1:10 %se definen 10 pasos para que se vea animada
    for e = 1:nelem

        m = elem(e, X); % Nodo i
        n = elem(e, Y); % Nodo j
        
        x_i = coord(m, X); % Extraigo la coordenada en X del nodo i
        y_i = coord(m, Y); % Extraigo la coordenada en Y del nodo i 
        x_j = coord(n, X); % Extraigo la coordenada en X del nodo j
        y_j = coord(n, Y); % Extraigo la coordenada en Y del nodo j
        
        xd_i = x_i + D(2*m-1)/10*paso; % Le sumo a la coordenada X del nodo i su desplazamiento (osea deformada)
        yd_i = y_i + D(2*m)/10*paso; % Le sumo a la coordenada Y del nodo i su desplazamiento (osea deformada)
        xd_j = x_j + D(2*n-1)/10*paso; % Le sumo a la coordenada X del nodo j su desplazamiento (osea deformada)
        yd_j = y_j + D(2*n)/10*paso; % Le sumo a la coordenada Y del nodo j su desplazamiento (osea deformada)
        
        h{e} = plot([xd_i, xd_j], [yd_i, yd_j], 'g');
    end
xlabel('x')
ylabel('y')
axis equal
title('CERCHA CON SU DEFORMADA')
pause(0.2) %detiene por 0.2 segundos el ciclo mientras se cumple la siguiente condicion
if paso < 10 %si no se ha llegado al paso final elimine la grafica de deformada anterior
    for e = 1:nelem
    delete(h{e}) 
    end
end
end

figure % Se añade otro figure para separar la figura de la cercha deformada de la cercha con los miembros a traccion, compresion y fuerza cero
for e = 1:nelem
    hold on
    m = elem(e, X); % Nodo i
    n = elem(e, Y); % Nodo j

    x_i = coord(m, X); % Extraigo la coordenada en X del nodo i
    y_i = coord(m, Y); % Extraigo la coordenada en Y del nodo i 
    x_j = coord(n, X); % Extraigo la coordenada en X del nodo j
    y_j = coord(n, Y); % Extraigo la coordenada en Y del nodo j

    color = 'k'; % "k" es la denominacion para el color negro que se usa en fuerza cero y como color inicial
    if N(e) > 0 % Los elementos a traccion se dibujan en azul en la figura 2
        color= 'b'; % "b" es la denominacion para el color azul por "blue"
    elseif N(e) < 0 % Los elementos a traccion se dibujan en rojo en la figura 2
        color = 'r';% "r" es la denominacion para el color rojo por "red"
    end
    plot([x_i, x_j], [y_i, y_j], color) % Comando para dibujar lo que esta en el ciclo
end
xlabel('x')
ylabel('y')
axis equal
title('GRAFICO A COLOR DE TRACCION, COMPRESION Y FUERZA CERO')

%% ************************ MUESTRA DE RESULTADOS ************************

format shortG
disp('-------- RESULTADOS --------')
disp('')

% MOSTRAR REACCIONES
disp('Estas son las Reacciones en los apoyos')
disp(Pa)
disp('')

%MOSTRAR DESPLAZAMIENTOS
disp('Estos son los desplazamientos de cada nodo')
disp(D)
disp('')

%MOSTRAR FUERZAS AXIALES DE CADA ELEMENTO
disp('Estas son las fuerzas axiales de cada elemento')
disp(N)