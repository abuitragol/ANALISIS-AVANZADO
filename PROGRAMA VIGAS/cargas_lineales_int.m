% f(x) = Ax^2+Bx+C

clc 
clear 
close all

 % Defininicion de variable simbolica

 syms x A B 
 
 %funcion
 
 carga_lineal = -A*x-B;

 % Integrar la funcion

 F1 = int(carga_lineal);
 F2 = int(F1);
 F3 = int(F2);
 F4 = int(F3);



% vector con las integrales de la funcion de la carga
F = [carga_lineal F1 F2 F3 F4]';
disp(F)


