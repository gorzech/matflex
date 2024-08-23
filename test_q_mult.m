% Some basic test for q_mult

p1 = normunit(rand(4, 1));
p2 = normunit(rand(4, 1));

pn = q_mult(p1, p2);

% Now the same with rotation matrices

A1 = Rot(p1);
A2 = Rot(p2);

An = Rot(pn);

A_expected = A1 * A2;

assert(norm(An - A_expected) < 1e-15)
