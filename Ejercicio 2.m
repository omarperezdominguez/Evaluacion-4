%Limpieza de pantalla
clear all
close all
clc

%SECCIÓN 1
%Declaración de variables simbólicas
syms th1(t) th2(t) th3(t) th4(t) th5(t) th6(t) t l0 l1 l2 l3 l4 l5 h

%SECCIÓN 2
%Configuración del robot, 0 para junta rotacional, 1 para junta prismática
RP = [0 0 0 0 0 0];

%SECCIÓN 3
%Creamos el vector de coordenadas articulares
Q = [th1(t), th2(t), th3(t), th4(t), th5(t), th6(t)];
disp('Coordenadas generalizadas');
pretty (Q);

%SECCIÓN 4
%Creamos el vector de velocidades generalizadas
Qp = diff(Q, t);
disp('Velocidades generalizadas');
pretty (Qp);

%SECCIÓN 5
%Número de grado de libertad del robot
GDL = size(RP, 2);
GDL_str = num2str(GDL);

%SECCIÓN 6
P = sym(zeros(3,1,GDL));
R = sym(zeros(3,3,GDL));




%Junta 1
%Posición de la junta 1 respecto a 0
P(:,:,1)= [0; 0; l0];
%Matriz de rotación de la junta 1 respecto a 0
Rz1 = [cos(th1)  0    sin(th1)  ;
       sin(th1)   0  -cos(th1);
       0       1       0];
R(:,:,1) = Rz1;



%Junta 2
%Posición de la junta 2 respecto a 1
P(:,:,2)= [-l1*sin(th2); l1*cos(th2); 0];
%Matriz de rotación de la junta 2 respecto a 1
Rz2 = [cos(th2) sin(th2)  0;
       sin(th2)  -cos(th2)  0;
       0         0         -1];
R(:,:,2) = Rz2;




%Junta 3
%Posición de la junta 3 respecto a 2
P(:,:,3) = [l2*sin(th3); l2*cos(th3); 0];
%Matriz de rotación de la junta 3 respecto a 2
Rz3 = [cos(th3) 0  sin(th3);
       sin(th3) 0 -cos(th3);
       0         1         0];
R(:,:,3) = Rz3;







%Junta 4
%Posición de la junta 4 respecto a 3
P(:,:,4)= [0; 0; l3];
%Matriz de rotación de la junta 4 respecto a 3
Rz4 = [cos(th4) 0 sin(th4) ;
       sin(th4)  0  -cos(th4);
      0         1       0];
R(:,:,4) = Rz4;








%Junta 5
%Posición de la junta 5 respecto a 4
P(:,:,5)= [-l4*sin(th5); l4*cos(th5); 0];
%Matriz de rotación de la junta 5 respecto a 4
Rz5 = [0 -sin(th5)  cos(th5);
       0  cos(th5)  sin(th5);
       -1         0         0];

R(:,:,5) = Rz5;









%Junta 6 
%Posición de la junta 6 respecto a 5
P(:,:,6) = [0; 0; 0];
%Matriz de rotación de la junta 6 respecto a 5
R(:,:,6) = [cos(th6) -sin(th6) 0;
            sin(th6)  cos(th6) 0;
            0         0        1];


%SECCIÓN 7
%Creamos un vector de ceros
Vector_Zeros= zeros(1, 3);

%Inicializamos las matrices de transformación Homogénea locales
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las matrices de transformación Homogénea globales
T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
%Inicializamos las posiciones vistas desde el marco de referencia inercial
PO(:,:,GDL)= P(:,:,GDL); 
%Inicializamos las matrices de rotación vistas desde el marco de referencia inercial
RO(:,:,GDL)= R(:,:,GDL); 
%Inicializamos las INVERSAS de las matrices de rotación vistas desde el marco de referencia inercial
RO_inv(:,:,GDL)= R(:,:,GDL); 

%SECCIÓN 8
for i = 1:GDL
    i_str= num2str(i);
    %Locales
    disp(strcat('Matriz de Transformación local A', i_str));
    A(:,:,i)=simplify([R(:,:,i) P(:,:,i); Vector_Zeros 1]);
    pretty (A(:,:,i));

    %Globales
    try
       T(:,:,i)= T(:,:,i-1)*A(:,:,i);
    catch
       T(:,:,i)= A(:,:,i);
    end
    disp(strcat('Matriz de Transformación global T', i_str));
    T(:,:,i)= simplify(T(:,:,i));
    pretty(T(:,:,i))
    
    RO(:,:,i)= T(1:3,1:3,i);
    RO_inv(:,:,i)= transpose(RO(:,:,i));
    PO(:,:,i)= T(1:3,4,i);
    %pretty(RO(:,:,i));
    %pretty(RO_inv(:,:,i));
    %pretty(PO(:,:,i));
end

%Calculamos la matriz de transformación del marco de referencia inercial
%visto desde el actuador final
% disp(strcat('Matriz de Transformación T', GDL_str,'_O calculada de forma manual'));
% RF_O=RO_inv(:,:,GDL);
% PF_O=-RF_O*PO(:,:,GDL);
% TF_O= simplify([RF_O PF_O; Vector_Zeros 1]);
% pretty(TF_O);   

%disp(strcat('Matriz de Transformación T', GDL_str,'_O calculada de forma auntomática'));
%pretty(simplify(inv(T(:,:,GDL))));

%SECCIÓN 9
%Calculamos el jacobiano lineal de forma diferencial
disp('Jacobiano lineal obtenido de forma diferencial');

% Inicializamos el Jacobiano
jv_d = sym(zeros(3, GDL));

% Calculamos derivadas parciales de x, y, z respecto a th1..th6
for i = 1:6
    jv_d(1,i) = functionalDerivative(PO(1,1,GDL), eval(['th', num2str(i)]));
    jv_d(2,i) = functionalDerivative(PO(2,1,GDL), eval(['th', num2str(i)]));
    jv_d(3,i) = functionalDerivative(PO(3,1,GDL), eval(['th', num2str(i)]));
end

% Simplificamos y mostramos
jv_d = simplify(jv_d);
pretty(jv_d);

%SECCIÓN 10
%Calculamos el jacobiano lineal de forma analítica
Jv_a(:,GDL) = PO(:,:,GDL);
Jw_a(:,GDL) = PO(:,:,GDL);

for k = 1:GDL
    if RP(k) == 0 %Casos: articulación rotacional
       %Para las juntas de revolución
        try
            Jv_a(:,k)= cross(RO(:,3,k-1), PO(:,:,GDL)-PO(:,:,k-1));
            Jw_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)= cross([0,0,1], PO(:,:,GDL));
            Jw_a(:,k)=[0,0,1];
        end
        
        %Para las juntas prismáticas
     elseif RP(k) == 1 %Casos: articulación prismática
%         
        try
            Jv_a(:,k)= RO(:,3,k-1);
        catch
            Jv_a(:,k)=[0,0,1];
        end
            Jw_a(:,k)=[0,0,0];
     end
 end    

Jv_a = simplify (Jv_a);
Jw_a = simplify (Jw_a);
disp('Jacobiano lineal obtenido de forma analítica');
pretty (Jv_a);
disp('Jacobiano ángular obtenido de forma analítica');
pretty (Jw_a);

disp('Velocidad lineal obtenida mediante el Jacobiano lineal');
V = simplify (Jv_a*Qp');
pretty(V);
disp('Velocidad angular obtenida mediante el Jacobiano angular');
W = simplify (Jw_a*Qp');

pretty(W);