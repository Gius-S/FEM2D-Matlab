function zh = phih_Pk(k, i, xh, yh)

% -------------------------------------------------------------------------
% Funzioni di base sul triangolo di riferimento Eh 
%
% -------------------------------------------------------------------------
%
% INPUT:
% 1. k: grado polinomiale
% 2. i: funzione di base scelta i = 1, \dots, (k+1)*(k+2)/2
% 3. [xh, yh]: punto in cui valutare la funzione
%
% OUTPUT
% 1. zh: valutazione della funzione bi base i-esima in [xh, yh]
%
% -------------------------------------------------------------------------

switch k
    
    % ---------------------------------------------------------------------
    % P1
    % ---------------------------------------------------------------------
    case 1
    switch i
        case 1
            zh = 1 - xh - yh;
        case 2
            zh = xh;
        case 3
            zh = yh;
    end
    
    % ---------------------------------------------------------------------
    % P2
    % ---------------------------------------------------------------------
    case 2
    switch i
        case 1
            zh = (1 - xh - yh).*(1 - 2*xh - 2*yh);
        case 2
            zh = xh.*(2*xh - 1);
        case 3
            zh = yh.*(2*yh - 1);
        case 4
            zh = 4*xh.*yh;
        case 5
            zh = -4*yh.*(xh + yh - 1);
        case 6
            zh = -4*xh.*(xh + yh - 1);
    end
    
    % ---------------------------------------------------------------------
    % P3
    % ---------------------------------------------------------------------
    case 3
    switch i
        case 1
            zh = 0.5*(1 - xh - yh).*(2 - 3*xh - 3*yh).*(1 - 3*xh - 3*yh);
        case 2
            zh = 0.5*xh.*(2 - 3*xh).*(1 - 3*xh);
        case 3
            zh = 0.5*yh.*(2 - 3*yh).*(1 - 3*yh);
        case 4
            zh = -9/2*xh.*yh.*(1 - 3*xh);
        case 5
            zh = -9/2*yh.*(1 - xh - yh).*(1 - 3*yh);
        case 6
            zh = 9/2*xh.*(1 - xh - yh).*(2 - 3*xh - 3*yh);
        case 7
            zh = -9/2*xh.*yh.*(1 - 3*yh);
        case 8
            zh = 9/2*yh.*(1 - xh - yh).*(2 - 3*xh - 3*yh);
        case 9
            zh = -9/2*xh.*(1 - xh - yh).*(1 - 3*xh);
        case 10
            zh = 27.*xh.*yh.*(1 - xh - yh);

    end
    
end

end