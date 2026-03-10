function [u_h, err] = FEM_2D(mesh, f_fun, k_fun, u_ex, du_ex)

% -------------------------------------------------------------------------
% Risolve il problema di diffusione stazionario
%
%    -div( kappa * grad u ) = f    in Omega
%                         u = u_ex su d_Omega
%
% -------------------------------------------------------------------------
% INPUT:
%   mesh    - struct con campi: xv, yv, vertices, edges, endpoints,
%             boundary, boundedges  (formato meshtrans.m)
%   f_fun   - function handle  f(x,y)         termine noto
%   k_fun   - function handle  kappa(x,y)     coefficiente di diffusione
%   u_ex    - function handle  u(x,y)         soluzione esatta
%   du_ex   - function handle  [dxu,dyu] = du_ex(x,y)  gradiente esatto
%
% OUTPUT:
%   u_h     - vettore (NV x 1) con i valori della soluzione FEM P1
%             ai vertici della mesh (inclusi i nodi di bordo)
%   err     - struct con campi:
%               err.L2    errore assoluto in norma L^2
%               err.H1    errore assoluto in seminorma H^1
%               err.Linf  errore assoluto in norma L^infty (ai nodi)
%

% =========================================================================
% 1. info dello struct mesh
% =========================================================================
xv         = mesh.xv;
yv         = mesh.yv;
vertices   = mesh.vertices;
boundary   = mesh.boundary;

NV  = length(xv);           % numero di vertici
NT  = size(vertices, 1);    % numero di triangoli
pik = 3;                    % dim(P1) in 2D: (1+1)*(1+2)/2 = 3

[xhq_A, yhq_A, whq_A] = quadratura(2);
Nq_A = length(whq_A);

[xhq_E, yhq_E, whq_E] = quadratura(4);
Nq_E = length(whq_E);


phih_A  = zeros(Nq_A, pik);
gxph_A  = zeros(Nq_A, pik);
gyph_A  = zeros(Nq_A, pik);
for i = 1:pik
    phih_A(:, i)               = phih_Pk(1, i, xhq_A, yhq_A);
    [gxph_A(:,i), gyph_A(:,i)] = grad_phih_Pk(1, i, xhq_A, yhq_A);
end


phih_E  = zeros(Nq_E, pik);
gxph_E  = zeros(Nq_E, pik);
gyph_E  = zeros(Nq_E, pik);
for i = 1:pik
    phih_E(:, i)               = phih_Pk(1, i, xhq_E, yhq_E);
    [gxph_E(:,i), gyph_E(:,i)] = grad_phih_Pk(1, i, xhq_E, yhq_E);
end


A = sparse(NV, NV);
F = zeros(NV, 1);

for ell = 1:NT

    v1 = vertices(ell, 1);
    v2 = vertices(ell, 2);
    v3 = vertices(ell, 3);

    x1 = xv(v1);  y1 = yv(v1);
    x2 = xv(v2);  y2 = yv(v2);
    x3 = xv(v3);  y3 = yv(v3);

    % Jacobiana della mappa affine
    J    = [x2-x1, x3-x1;
        y2-y1, y3-y1];
    JIT  = inv(J)';
    area = 0.5 * abs(det(J));
    xB    = (x1 + x2 + x3) / 3;
    yB    = (y1 + y2 + y3) / 3;
    k_ell = k_fun(xB, yB);

    AE = zeros(pik, pik);
    FE = zeros(pik, 1);

    for i = 1:pik
        for j = i:pik
            for q = 1:Nq_A
                grad_j = JIT * [gxph_A(q, j); gyph_A(q, j)];
                grad_i = JIT * [gxph_A(q, i); gyph_A(q, i)];
                AE(i, j) = AE(i, j) + dot(grad_j, grad_i) * whq_A(q);
            end
            AE(i, j) = 2 * area * k_ell * AE(i, j);
            AE(j, i) = AE(i, j);   % simmetria
        end


        for q = 1:Nq_A
            pq = J * [xhq_A(q); yhq_A(q)] + [x1; y1];
            FE(i) = FE(i) + f_fun(pq(1), pq(2)) * phih_A(q, i) * whq_A(q);
        end
        FE(i) = 2 * area * FE(i);

    end

    INDICI = [v1; v2; v3];
    A(INDICI, INDICI) = A(INDICI, INDICI) + AE;
    F(INDICI)         = F(INDICI) + FE;

end


u_h = zeros(NV, 1);

u_h(boundary) = u_ex(xv(boundary), yv(boundary));

internV = setdiff(1:NV, boundary);

F_int = F(internV) - A(internV, boundary) * u_h(boundary);

A_int = A(internV, internV);

u_h(internV) = A_int \ F_int;

err.L2   = 0;
err.H1   = 0;
err.Linf = 0;

uex_nodi = u_ex(xv, yv);
err.Linf = max(abs(uex_nodi - u_h));

errL2_sq = 0;
errH1_sq = 0;

for ell = 1:NT

    v1 = vertices(ell, 1);
    v2 = vertices(ell, 2);
    v3 = vertices(ell, 3);

    x1 = xv(v1);  y1 = yv(v1);
    x2 = xv(v2);  y2 = yv(v2);
    x3 = xv(v3);  y3 = yv(v3);

    J    = [x2-x1, x3-x1;
        y2-y1, y3-y1];
    JIT  = inv(J)';
    area = 0.5 * abs(det(J));

    INDICI = [v1; v2; v3];
    uh_ell = u_h(INDICI);

    for q = 1:Nq_E
        pq = J * [xhq_E(q); yhq_E(q)] + [x1; y1];
        xq = pq(1);
        yq = pq(2);
        uh_q  = 0;
        guh_q = [0; 0];
        for i = 1:pik
            uh_q  = uh_q  + uh_ell(i) * phih_E(q, i);
            guh_q = guh_q + uh_ell(i) * (JIT * [gxph_E(q,i); gyph_E(q,i)]);
        end

        uex_q          = u_ex(xq, yq);
        [dxu_q, dyu_q] = du_ex(xq, yq);
        guex_q         = [dxu_q; dyu_q];
        errL2_sq = errL2_sq + 2*area * (uex_q - uh_q)^2       * whq_E(q);
        errH1_sq = errH1_sq + 2*area * norm(guex_q - guh_q)^2 * whq_E(q);

    end

end

err.L2 = sqrt(errL2_sq);
err.H1 = sqrt(errH1_sq);

end