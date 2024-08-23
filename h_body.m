function body = h_body(mass, Ic)
%H_BODY Simple helper that returns new body

assert(mass > 0, "Mass must be positive");
if isvector(Ic)
    Ic = diag(Ic);
end
assert(all(size(Ic) == [3, 3]), "Wrong size of an inertia matrix")
assert(all(diag(Ic) > 0), "Inertia diagonal terms must be positive")

body.mass = mass;
body.Ic = Ic;

body.nh = 6;
body.nq = 7;
end

