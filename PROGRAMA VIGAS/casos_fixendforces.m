%% ------------- Tipos de carga -------------
% datos de entrada
W = 0; % Carga puntual
l_1 = 0; % Distancia hasta el nodo izquierdo
l_2 = 0; % Distancia hasta el nodo derecho
L = 0; % longitud total del elemento
M = 0; % momento de par
w = 0; % carga distribuida
%   |------l_1-----w ------l_2----|    (distancias de la carga a los nodos)
%                  |
%                  v 
%  b|-----------------------------|e   (viga con ndoos b y e y fuerza (W) )
%   |--------------L--------------|    (longitud total de la viga (L) )
%% Carga puntual
Fsb = ((W*l_2^2)/L^3)+(3*l_1+l_2); % reaccion en b
Fmb = ((W*l_1*l_2^2)/L^2); % momento con respecto a b
Fse = ((W*l_1^2)/L^3)+(l_1+3*l_2); % reaccion en e
Fme = -((W*l_1^2*l_2)/L^2); % momento con respecto a e

%% Momento de par
Fsb = -(6*M*l_1*l_2)/L^3; % reaccion en b
Fmb = ((M*l_2)/L^2)*(l_2-2*l_1); % momento con respecto a b
Fse = (6*M*l_1*l_2)/L^3; % reaccion en e
Fmb = ((M*l_1)/L^2)*(l_1-2*l_2); % momento con respecto a e

%% Carga distribuida 

Fsb = ((w*L)/2)*(1-(((l_1/L^4)*(2*L^3-2*l_1^2*L+l_1^3))-((l_2^3/L^4)*(2*L-L_2)))); % reaccion en b
Fmb = ((w*l^2)/12)*(1-((l_1^2/L^4)*(6*L^2-8*l_1*L+3*l_1^2))-((l_2^3/L^4)*(4*L-3*l_2))); % momento con respecto a b
Fse = ((w*L)/2)*(1-(((l_1^3)/L^4)*2*L-1)-((l_2^2/L^4)*(2*L^3-2*l_2^2*L+l_2^3))); % reaccion en e
Fme = -((w*l^2)/12)*(1-((l_1^3/L^4)*(4*L-3*l_1)-((l_2^2/L^4)*(6*L^2-8*l_2*L+3*l_2^2)))); % momento con respecto a e
