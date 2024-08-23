function A = Rot(p)
% A = Rot(p)
% funkcja pobiera parametry eulera (p) i podaje dla nich macierz rotacji A

A = Eep(p)*Gep(p)';