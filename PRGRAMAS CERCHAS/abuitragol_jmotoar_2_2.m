%FALTA ACOMODAR LAS FUERZAS EN EL APOYO DE LA DERECHA PARA QUE QUEDEN CON
%EL ANGULO DEL RODILLO INCLINADO
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
a = [1 2 12];

% Grados de libertad en el nodo con apoyo inclinado

a_incl = [11 12];

% angulo de inclinacion del apoyo en grados

alpha = 60;

% Despazamientos gdl restringidos
Da = [0 0 0 ]';

% Grados de libertad no restringidos
b = [3 4 5 6 7 8 9 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 ];

% Fuerzas en los gdl no restringidos (fuerzas nodales externas)
Pb = [0 0 0 0 0 0 0 0 0 0 -20*P -(20*P)*cosd(69.31) -(20*P)*sind(69.31) 0 -30*P 0 0 0 0 0 -30*P (40*P)*cosd(69.31) -(40*P)*sind(69.31) 0 -40*P 0 -50*P ]';

% Fuerzas en los gdl de los apoyos 
Fa = [0 -10*P -10*P]';
%% ------ (2) PROCESO DE CALCULO --------

% Se separa la memoria en celdas que son contenedores (cajones) que pueden
% contener cualquier tipo de dato

ke = cell(nelem,1);

T = cell(nelem,1);

Ke = cell(nelem,1);

K = zeros(ngdl,ngdl);

% ****** ENSAMBLAJE MATRICIAL ******

% Con ciclo for para optimizar el proceso ya que es un proceso repetitivo
% Matriz de rigidez ke en sistema local Eq(17)

for e = 1:nelem
ke = [A*E/L(e) 0 -A*E/L(e) 0;0 0 0 0; -A*E/L(e) 0 A*E/L(e) 0; 0 0 0 0];

% MatriZ de transformacion de cada elemento Eq(22)
T = [eta(e) mu(e) 0 0; -mu(e) eta(e) 0 0; 0 0 eta(e) mu(e); 0 0 -mu(e) eta(e)];

% Matriz de rigidez Ke en sistema global Eq(25)
Ke = T'*ke*T;

% Grados de libertad globales del elemento e
gdl_e = LaG(e,:);

% Se suma las contribuciones ke a la matriz K
K(gdl_e, gdl_e) = K(gdl_e, gdl_e) + Ke;

end

%% ----------- Matriz de rigidez mixta debido al apoyo inclinado -----------

%Matriz de transformacion espacial
Tm = eye(ngdl,ngdl); % eye crea una matriz identidad de dimensiones dadas
Tm(a_incl,a_incl)= [cosd(alpha) sind(alpha); -sind(alpha) cosd(alpha)]; % sind y cosd ya que el angulo esta en grados "degrees"

%Matriz de rigidez  mixta
Km = Tm*K*Tm';

% Vectores de fuerza y desplazamiento mixto 

% Se generan los vectores D, P y F
D = zeros(ngdl,1);
D(a) = Da;

P = zeros(ngdl,1);
P(b) = Pb;

F = zeros(ngdl,1);
F(a) = Fa;

%se vuelven mixtos
Dm = Tm*D;
Pm = Tm*P;
Fm = Tm*F;

%Se extrae el Da mixto y el Pb mixto
Dam = Dm(a);
Pbm = Pm(b);
Fam = Fm(a);


% Particion de las submatrices en el sistema mixto
% ( Pa )     ( Kaa | Kab )  (Da)
%-------  =  ------|------  ----
% ( Pb )     ( Kba | Kbb )  (Db)

Kaa = Km(a,a);
Kab = Km(a,b);
Kba = Km(b,a);
Kbb = Km(b,b);


%se calculan los desplazamientos desconocidos Eq(7) y Eq(8)

Dbm = Kbb\(Pb-Kba*Dam);
Pam = (Kaa*Dam + Kab*Dbm) - Fam;
% Y se guardan los vectores Dm y Pm

    Dm(b) = Dbm; % Desplazamientos en el sistema mixto
    Pm(a) = Pam; % Reacciones en el sistema mixto 
    Fm(a) = Fam;
% Se regresan los resultados al sistema global D y P

D = Tm'*Dm;
Da = D(union(a,a_incl));
Db = D(b);

P = Tm'*Pm;
Pa = P(union(a,a_incl));
Pb = P(b);

F = Tm'*Fm;
Fa = F(union(a,a_incl));
Fa = F(a);
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