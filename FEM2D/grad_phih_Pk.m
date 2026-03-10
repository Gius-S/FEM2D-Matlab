function [gxh, gyh] = grad_phih_Pk(k, i, xh, yh)

% -------------------------------------------------------------------------
% Gradiente delle funzioni di base sul triangolo di riferimento Eh 
%
% -------------------------------------------------------------------------
%
% INPUT:
% 1. k: grado polinomiale
% 2. i: funzione di base scelta i = 1, \dots, (k+1)*(k+2)/2
% 3. [xh, yh]: punto in cui valutare la funzione
%
% OUTPUT
% 1. gxh: valutazione della componente x della funzione di base i-esima in [xh, yh]
% 2. gxy: valutazione della componente y della funzione di base i-esima in [xh, yh]
% -------------------------------------------------------------------------

switch k
    
    % ---------------------------------------------------------------------
    % P1
    % ---------------------------------------------------------------------
    case 1
    switch i
        case 1
            gxh = -ones(size(xh));
            gyh = -ones(size(xh));
        case 2
            gxh = ones(size(xh));
            gyh = zeros(size(xh));
        case 3
            gxh = zeros(size(xh));
            gyh = ones(size(xh));
    end
    
    % ---------------------------------------------------------------------
    % P2
    % ---------------------------------------------------------------------
    case 2
    switch i
        case 1
            gxh = 4*xh + 4*yh - 3;
            gyh = 4*xh + 4*yh - 3;
        case 2
            gxh = 4*xh - 1;
            gyh = zeros(size(xh));
        case 3
            gxh = zeros(size(xh));
            gyh = 4*yh - 1;
        case 4
            gxh = 4*yh;
            gyh = 4*xh;
        case 5
            gxh = -4*yh;
            gyh =  4 - 8*yh - 4*xh;
        case 6
            gxh = 4 - 4*yh - 8*xh;
            gyh = -4*xh;
    end
    
    % ---------------------------------------------------------------------
    % P3
    % ---------------------------------------------------------------------
    case 3
    switch i
        case 1
            gxh = 18*xh + 18*yh - 27*xh.*yh - (27*xh.^2)/2 - (27*yh.^2)/2 - 11/2;
            gyh = 18*xh + 18*yh - 27*xh.*yh - (27*xh.^2)/2 - (27*yh.^2)/2 - 11/2;
        case 2
            gxh = (27*xh.^2)/2 - 9*xh + 1;
            gyh = zeros(size(xh));
        case 3
            gxh = zeros(size(xh));
            gyh = (27*yh.^2)/2 - 9*yh + 1;
        case 4
            gxh = (9*yh.*(6*xh - 1))/2;
            gyh = (9*xh.*(3*xh - 1))/2;
        case 5
            gxh = -(9*yh.*(3*yh - 1))/2;
            gyh = (9*xh)/2 + 36*yh - 27*xh.*yh - (81*yh.^2)/2 - 9/2;
        case 6
            gxh = (81*xh.^2)/2 + 54*xh.*yh - 45*xh + (27*yh.^2)/2 - (45*yh)/2 + 9;
            gyh = (9*xh.*(6*xh + 6*yh - 5))/2;
        case 7
            gxh = (9*yh.*(3*yh - 1))/2;
            gyh = (9*xh.*(6*yh - 1))/2;
        case 8
            gxh = (9*yh.*(6*xh + 6*yh - 5))/2;
            gyh = (27*xh.^2)/2 + 54*xh.*yh - (45*xh)/2 + (81*yh.^2)/2 - 45*yh + 9;
        case 9
            gxh = 36*xh + (9*yh)/2 - 27*xh.*yh - (81*xh.^2)/2 - 9/2;
            gyh = -(9*xh.*(3*xh - 1))/2;
        case 10
            gxh = -27*yh.*(2*xh + yh - 1);
            gyh = -27*xh.*(xh + 2*yh - 1);
    end
    
end

end