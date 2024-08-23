function est = sksym(e)
% funkcja liczaca macierz skosnosymetryczna do wektora e(3x1)
est = [0 -e(3) e(2); 
    e(3) 0 -e(1); 
    -e(2) e(1) 0];