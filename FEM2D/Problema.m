% Script Matlab per la risoluzione del problema:  - Δu = f    con 
%
%           u = |r|^2/3 sin(2/3 theta)      
%           f = 0 
%
%   Nel dominio L = (-1,1)x(-1,1) \ (-1,0]x(-1,0]

clear; clc; close all;

%%  DEFINIZIONE DEL PROBLEMA 

syms x y real

r = sqrt(x^2 + y^2);
theta = atan2(y, x);

u_polare = r^(2/3) * sin((2*theta+ pi)/3);
u_sym= simplify(u_polare); % funzione 

dux_sym = diff(u_sym, x);
duy_sym = diff(u_sym, y);


u_ex = matlabFunction(u_sym, 'Vars', [x, y]);
f_fun = @(x,y) zeros(size(x));   % termine noto nullo
dux_fun = matlabFunction(dux_sym, 'Vars', [x, y]);
duy_fun = matlabFunction(duy_sym, 'Vars', [x, y]);
du_ex = @(x,y) deal(dux_fun(x,y), duy_fun(x,y));

k_fun = @(x,y) ones(size(x)); % coefficente di diffusione 


%% Inizializzazione variabili 

% Studio Errore 
h_vals = logspace(log10(0.2), log10(0.025), 20);  
NCasi = length(h_vals);

errorL2 = zeros(1, NCasi);
errorH1 = zeros(1, NCasi);

h_eff = zeros(1, NCasi);

%Dominio

R1 = [3; 4; -1; 1; 1; -1;  0; 0; 1; 1];   
R2 = [3; 4;  0; 1; 1;  0; -1; -1; 0; 0]; 
Q = [R1, R2];
sf = 'R1+R2';
ns = char('R1', 'R2')';
D = decsg(Q, sf, ns);


%% Approssimazione 

for i = 1 : NCasi
    fprintf('Solving for h = %f ... ', h_vals(i));

    % Mesh 
    [P, E, T] = initmesh(D, 'Hmax', h_vals(i));
    mesh = meshtrans(P, E, T);
    area_tot = 0; 
    for e = 1:size(mesh.vertices, 1) %itero sui triangoli
        v = mesh.vertices(e, :); %info sui vertici 
        x1 = mesh.xv(v(1)); y1 = mesh.yv(v(1)); 
        x2 = mesh.xv(v(2)); y2 = mesh.yv(v(2));
        x3 = mesh.xv(v(3)); y3 = mesh.yv(v(3));
        area_tot = area_tot + abs(0.5 * det([x2-x1, x3-x1; y2-y1, y3-y1])); %calcolo area triangolo
    end
    h_eff(i) = sqrt(area_tot / size(mesh.vertices, 1));

    [u_num, err] = FEM_2D(mesh, f_fun, k_fun, u_ex, du_ex); % Soluzione Numerica 
    
    errorL2(i) = err.L2;
    errorH1(i) = err.H1;

    fprintf('Done.\n');
end

%% Print Risultati 


figure(1);

loglog(h_eff, errorL2, 'bo-', 'LineWidth', 2, 'MarkerSize', 7);
hold on;
loglog(h_eff, errorH1, 'gs-', 'LineWidth', 2, 'MarkerSize', 8);

% Linee di riferimento teoriche
loglog(h_eff, h_eff.^(4/3), 'y--', 'LineWidth', 1.5);   % O(h^{4/3}) per L2
loglog(h_eff, h_eff.^(2/3), 'y:', 'LineWidth', 1.5);    % O(h^{2/3}) per H1

xlabel('h (passo medio)');
ylabel('Errore');
legend('L^2', 'H^1', 'O(h^{4/3})', 'O(h^{2/3})', 'Location', 'southeast');
title('Convergenza FEM P1 - Dominio a L');
xlim([1.1763e-02, 1.3043e-01]);
ylim([2.9317e-04, 5.4666e-01]);
zlim([-1.0000e+00, 1.0000e+00]);

pL2 = polyfit(log(h_eff), log(errorL2), 1);

pH1 = polyfit(log(h_eff), log(errorH1), 1);


fprintf('\n======== STIMA ORDINI ========\n\n');
fprintf(['Ordine L2:\n' ...
         '              %.3f     expected: %.3f\n'], pL2(1), 4/3);
fprintf(['Ordine H1:\n ' ...
         '              %.3f     expected: %.3f\n'], pH1(1), 2/3);


figure(2);

trisurf(mesh.vertices, mesh.xv, mesh.yv, u_num);

xlabel('x'); ylabel('y'); zlabel('u');



title('Soluzione numerica FEM (dominio a L)');

colorbar;
view([323.148 24.070]);