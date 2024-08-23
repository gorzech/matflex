function pn = q_mult(p1, p2) 
%Q_MULT Multiplication of quaternions. For rotational quat. rotation p2 
%followed by p1: A(p1) * A(p2)
pn = [p1(1)*p2(1) - p1(2)*p2(2) - p1(3)*p2(3) - p1(4)*p2(4)
     p1(1)*p2(2) + p1(2)*p2(1) + p1(3)*p2(4) - p1(4)*p2(3)
     p1(1)*p2(3) - p1(2)*p2(4) + p1(3)*p2(1) + p1(4)*p2(2)
     p1(1)*p2(4) + p1(2)*p2(3) - p1(3)*p2(2) + p1(4)*p2(1)];
 