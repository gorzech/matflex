function sys = sys_generate(bodies, fbodies, varargin)
%SYS_GENERATE Create sys structure as a single entry to a whole multibody
%system - needs to add more parts here later on

sys.nbodies = length(bodies);
sys.nfbodies = length(fbodies);
sys.fbodies = fbodies;
sys.bodies = bodies;

p = inputParser;
addParameter(p, 'grav', [0; 0; 0]);
addParameter(p, 'usequadratic', true);
% Baumgartes parameters
addParameter(p, 'alfa', 20);
addParameter(p, 'beta', 20);
parse(p, varargin{:});

sys.usequadratic = p.Results.usequadratic;
sys.baumg_p1 = 2 * p.Results.alfa;
sys.baumg_p2 = p.Results.beta ^ 2;

sys.grav = p.Results.grav;
sys.isgrav = norm(sys.grav) > 0;

sys.nconstr = 0;
sys.constr = struct('simple', [], 'point', [], 'perp1', [], ...
    'drive_point', []);
sys.force = struct('point', [], 'tsd', [], 'torque', []);
sys.nh = 6 * sys.nbodies;
sys.nq = 7 * sys.nbodies;

h_id = 0;
q_id = 0;
for ii = 1 : sys.nbodies
    sys.bodies(ii).h_idx = h_id + (1 : 6);
    sys.bodies(ii).q_idx = q_id + (1 : 7);
    h_id = h_id + 6;
    q_id = q_id + 7;
end

for ii = 1 : sys.nfbodies
    sys.nh = sys.nh + fbodies(ii).nh;
    sys.nq = sys.nq + fbodies(ii).nq;
    
    sys.fbodies(ii).h_idx = h_id + (1 : fbodies(ii).nh);
    sys.fbodies(ii).q_idx = q_id + (1 : fbodies(ii).nq);
    h_id = h_id + fbodies(ii).nh;
    q_id = q_id + fbodies(ii).nq;
    
    % What to run - depends on the flexible body type.
    if isfield(sys.fbodies(ii), 'm2')
        % Preintegrated ffr version
        sys.fbodies(ii).ffr_mass = @ffr_mass_preintegrated;
        sys.fbodies(ii).ffr_quad_inertia = @ffr_quadratic_inertia_preintegrated;
        sys.fbodies(ii).ffr_grav = @ffr_grav_force_preint;
    elseif isfield(sys.fbodies(ii), 'sec')
        % Version for integration
        sys.fbodies(ii).ffr_mass = @ffr_mass_integrate;
        sys.fbodies(ii).ffr_quad_inertia = @ffr_quadratic_inertia_simple_integrate;
        sys.fbodies(ii).ffr_grav = @ffr_grav_force;
    elseif isfield(sys.fbodies(ii), 'ksigma') 
        % SID body
        sys.fbodies(ii).ffr_mass = @ffr_mass_sid;
        sys.fbodies(ii).ffr_quad_inertia = @ffr_quadratic_inertia_sid;
        sys.fbodies(ii).ffr_grav = @ffr_grav_force_sid;
    else
        error('Unknown flexible body type at id %d', ii);
    end
end

sys.ny = sys.nh + sys.nq;
sys = rigid_body_mass(sys);
    
end

function sys = rigid_body_mass(sys)
bodies = sys.bodies;
    
i_body = [1, 2, 3, 4, 5, 6, 4, 5, 6, 4, 5, 6];
j_body = [1, 2, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6];
n_val = 12;

i_idx = zeros(1, n_val * sys.nbodies);
j_idx = zeros(size(i_idx));
r_val = zeros(size(i_idx));

b_idx = 0;
for ii = 0 : sys.nbodies - 1
    i_idx(b_idx + (1 : n_val)) = 6 * ii + i_body;
    j_idx(b_idx + (1 : n_val)) = 6 * ii + j_body;
    r_val(b_idx + (1 : 3)) = bodies(ii + 1).mass;
    r_val(b_idx + (4 : n_val)) = reshape(bodies(ii + 1).Ic, [1, 9]);
    b_idx = b_idx + n_val;
end

sys.m_rb = sparse(i_idx, j_idx, r_val);
end

