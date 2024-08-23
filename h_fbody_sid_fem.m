function fbody = h_fbody_sid_fem(file_name)

[sft, fbody] = init_read(file_name);

name = get_data(sft);
assert(strcmp(name, "new modal"), "Main class is called modal")

while true
    sft = advance(sft);
    name = get_data(sft);
    switch name
        case "end modal"
            break; % whole model read. Do not expect more
        case "refmod"
            [sft, fbody] = read_refmod(sft, fbody);
        case "frame"
            [sft, fbody] = read_frame(sft, fbody);
        otherwise
            [sft, fbody.(name)] = read_taylor(sft, name);
    end
    
end

sft = advance(sft);
name = get_data(sft);
assert(strcmp(name, "end part"), "SID file should end with part close tag")

fbody.nflex = fbody.n_mode;
fbody.nh = 6 + fbody.nflex;
fbody.nq = fbody.nh + 1;

% Need to rename Ke to Kff
fbody.Kff = fbody.Ke;
fbody = rmfield(fbody, 'Ke');

end

function [sft, fbody] = read_frame(sft, fbody)
while true
    sft = advance(sft);
    [name, value] = get_data(sft);
    switch name
        case "end frame"
            break;
        case "new node"
            [sft, fbody] = read_node(sft, fbody, value);
    end
end
end

function [sft, fbody] = read_node(sft, fbody, node_id)
while true
    sft = advance(sft);
    name = get_data(sft);
    switch name
        case "end node"
            break;
        case "origin"
            [sft, fbody.node(node_id).origin] = read_taylor(sft, name);
        case "phi"
            [sft, fbody.node(node_id).phi] = read_taylor(sft, name);
        case "psi"
            [sft, fbody.node(node_id).psi] = read_taylor(sft, name);
        case "AP"
            [sft, fbody.node(node_id).AP] = read_taylor(sft, name);
    end
end
end

function [sft, m] = read_taylor(sft, name)
% First approach - read all matrices into dense formats
end_name = sprintf('end %s', name);
order = 0;
nrow = 1;
ncol = 1;
nq = 0;
% nqn = 0;
structure = 0;
any_data_set = false;
while true
    sft = advance(sft);
    [name, value] = get_data(sft);
    switch name
        case end_name
            break;
        case "order"
            order = value;
            assert(order < 2, "Only Taylor matrices up to order 1 are supported")
        case "nrow"
            nrow = value;
        case "ncol"
            ncol = value;
        case "nq"
            nq = value;
        case "nqn"
            nqn = value;
            assert(nqn == 0, "Only Taylor matrices up to order 1 are supported")
        case "structure"
            structure = value;
        otherwise
            % Here only terms m0 and m1 are allowed
            assert(startsWith(name, 'm'), "Only matrix entries are allowed")
            if ~any_data_set
                m_tmp = init_data();
                any_data_set = true;
            end
            m_tmp = set_data(m_tmp, name, value, structure);
    end
end
if any_data_set
    if m_tmp.set_m0 
        m.m0 = m_tmp.m0;
    else
        m.m0 = [];
    end
    if m_tmp.set_m1
        m.m1 = m_tmp.m1;
    else
        m = m.m0;
    end
elseif structure == 4 % Unit matrix
    m = eye(nrow);
else    
    % Or m = []?
    m = zeros(nrow, ncol);
end
    function m_tmp = set_data(m_tmp, name, value, structure)
        if startsWith(name, 'm0')
            idx = sscanf(name, 'm0(%d,%d)', 2);
            m_tmp.set_m0 = true;
            m_tmp.m0(idx(1), idx(2)) = value;
            if structure == 2 % symmetric
                m_tmp.m0(idx(2), idx(1)) = value;
            end
        elseif startsWith(name, 'm1')
            idx = sscanf(name, 'm1(%d,%d,%d)', 3);
            m_tmp.set_m1 = true;
            m_tmp.m1(idx(1), idx(2), idx(3)) = value;
            if structure == 2 % symmetric
                m_tmp.m1(idx(3), idx(2), idx(1)) = value;
            end
        else
            error('Unknown matrix data entry')
        end
    end
    function m_tmp = init_data
        m_tmp.m0 = zeros(nrow, ncol);
        if order > 0
            m_tmp.m1 = zeros(nrow, nq, ncol);
        end
        m_tmp.set_m0 = false;
        m_tmp.set_m1 = false;
    end
end

function [sft, fbody] = read_refmod(sft, fbody)
while true
    sft = advance(sft);
    [name, value] = get_data(sft);
    switch name
        case "end refmod"
            break;
        case "mass"
            fbody.mass = value;
        case "nelastq"
            assert(value == fbody.n_mode, "Number of modes does not agree with refmod.nelastq")
    end    
end
end

function [sft, fbody] = init_read(file_name)
sft.data = sid_fem_raw_read(file_name);
sft.n_data = height(sft.data);
assert(sft.n_data > 3)

[A, n] = sscanf(sft.data.name(1), '%d', 2);

assert(n == 2, "Header should have 2 integer values")
assert(all(A > 0), "Counts in header should be positive integers")

fbody = struct('n_node', A(1), 'n_mode', A(2));
sft.idx = 2;

name = get_data(sft);
assert(strcmp(name, "part"), "File should begin with part definition")
name = get_data(sft, sft.n_data);
assert(strcmp(name, "end part"), "File should end with end part statement")
sft = advance(sft);
end

function sft = advance(sft)
if sft.idx + 1 > sft.n_data
    error('Cannot advance file read out of bounds.');
end
sft.idx = sft.idx + 1;
end

function [name, value] = get_data(sft, idx)
if nargin < 2
    idx = sft.idx;
end
if idx > sft.n_data
    error('Index out of bound in call to get_data.');
end
name = strtrim(sft.data.name(idx));
value = sft.data.value(idx);
end

function sid_fem_tbl = sid_fem_raw_read(filename)
%IMPORTFILE Import data from a text file
% Auto-generated by MATLAB on 18-Feb-2020 15:37:32
% Modified by Grzegorz Orzechowski on 20.02.2020

opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "=";

% Specify column names and types
opts.VariableNames = ["name", "value"];
opts.VariableTypes = ["string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "name", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "name", "EmptyFieldRule", "auto");

% Import the data
sid_fem_tbl = readtable(filename, opts);
end