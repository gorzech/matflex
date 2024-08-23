% Test if data generated for straight beam are as expected

%% One element case

n_elem = 1;
elem_length = 0.5;

[q, elem_idx, node_idx] = straight_beam3d(elem_length, n_elem);

assert(isvector(q))
assert(length(q) == 12)

q_expected = [0; 0; 0; 0; 0; 0
    elem_length; 0; 0; 0; 0; 0];
assert(norm(q_expected - q) < 1e-16)

assert(isvector(elem_idx))
assert(all(size(elem_idx) == [12, 1]))
assert(all(elem_idx == (1 : 12)'))

assert(all(size(node_idx) == [6, 2]))
assert(all(node_idx(:) == (1 : 12)'))

%% 4th element 

element_len = pi / 3;

q_4_expected = [element_len * 3; zeros(5, 1); element_len * 4; zeros(5, 1)];
ei_4_expected = (19:30)';
ni_45_expected = [19 : 24; 25 : 30]';

for ii = 4:8
    [q, elem_idx, node_idx] = straight_beam3d(element_len, ii);
    assert( all(size(q) == [6 * ii + 6, 1]), 'n_elem = %d\n', ii )
    assert( all(size(elem_idx) == [12, ii]), 'n_elem = %d\n', ii )
    assert( all(size(node_idx) == [6, ii + 1]), 'n_elem = %d\n', ii )
    assert( all(ei_4_expected == elem_idx(:, 4)), 'n_elem = %d\n', ii )
    assert( norm(q_4_expected - q(ei_4_expected)) < 1e-15, 'n_elem = %d\n', ii )
    assert( all(ni_45_expected == node_idx(:, [4, 5]), 'all'), 'n_elem = %d\n', ii )
end

