%% Test for ffr_grav_force - check if gravity vector is as expected

rng(31)

g = [1; 2; -9.81];

p = rand(4, 1);
p = p ./ norm(p);

file_name = 'test_import_straight_beam.SID_FEM';
fbody = h_fbody_sid_fem(file_name);

qf = rand(fbody.nflex, 1) * 0.173174;

M_sid = ffr_mass_sid(fbody, p, qf);

Qg = ffr_grav_force_sid(fbody, g, p, qf);

acc = M_sid \ Qg;

acc_expected = [g; zeros(3 + fbody.n_mode, 1)];

diff_acc = acc_expected - acc;
assert(norm(diff_acc) < 5e-14)

