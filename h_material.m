function mat = h_material(rho, E, nu)
%H_MATERIAL Helper function to create structure with material properties

if nargin < 1
    rho = 7801;
end
if nargin < 2
    E = 2e11;
end
if nargin < 3
    nu = 0.3;
end


mat.rho = rho;
mat.E = E;
mat.nu = nu;
mat.G = E / 2 / (1 + nu);

end

