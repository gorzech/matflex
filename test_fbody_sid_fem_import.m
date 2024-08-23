% Script to test if import of the SID_FEM seems to work properly

file_name = 'test_import_straight_beam.SID_FEM';

try
    fbody = h_fbody_sid_fem(file_name);
catch
    assert(false, 'Error when reading SID_FEM test file "%s"', file_name)
end

assert(isstruct(fbody) && isscalar(fbody), "fbody should be scalar struct")

assert(fbody.n_node == 2 && fbody.n_mode == 8)

assert(abs(fbody.mass - 9.585869E+01) < 1e-5)

assert(length(fbody.node) == 2)

assert(abs(fbody.node(2).AP.m1(1, 8, 3) - 2.26780545078E+00) < 1e-7)

assert(norm(fbody.Ct) == 0)

assert(isstruct(fbody.Oe))

assert(isnumeric(fbody.Kff))

assert(isequal(size(fbody.Kff), [8, 8]))