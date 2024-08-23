% Test if block overlap works as expected

% One block only - if this not works, why to go further?!

B_single = [2 0.5; 0 3];
B_idx = [1; 2];

B = block_overlap(B_idx, B_single);

assert( norm(B_single - B) == 0 )

%% Just two blocks 

B_idx = [1, 2; 2, 3];

B = block_overlap(B_idx, B_single);

B_expected = [2, 0.5, 0; 0, 5, 0.5; 0, 0, 3];

assert( norm(B_expected - B) < 1e-16 )

%% Few more blocks of random matrices

n_size = 5;
B_single = rand(5);
B_single(rand(5) > 0.5) = 0; % set some elements to zero

n_blocks = 6;

B_idx = [1:5
    4:2:12
    3:7
    11:15
    8:12
    6:10]';

% construct matrix by hand
B_expected = zeros(max(max(B_idx)));
for ii = 1 : n_blocks
    B_expected(B_idx(:, ii), B_idx(:, ii)) = ...
        B_expected(B_idx(:, ii), B_idx(:, ii)) + B_single;
end

B = block_overlap(B_idx, B_single);

assert( norm(B_expected - B) < 1e-16 )