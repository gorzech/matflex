function sec = h_section(parameters, type)
%H_SECTION Helper function that returns section properties - default is
%rectangular

if nargin < 1
    parameters = 0.05;
end
if nargin < 2
    type = 'rect';
end

switch type
    case 'rect'
        a = parameters(1);
        sec.A = a * a;
        sec.Iy = a ^ 4 / 12;
        sec.Iz = sec.Iy;
        sec.Jx = a ^ 4 / 6; % polar moment of inertia
        sec.thky = a;
        sec.thkz = a;
    otherwise
        error("Unknown section type %s!", type)
end

end

