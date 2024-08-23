function M = u_tay_mult(m, q)
%TAY_MULT General multiplication subroutine for Taylor matrices generated
%during SID_FEM import

if isnumeric(m) 
    M = m;
elseif isempty(m.m0)
    m1_mult();
elseif isempty(m.m1)
    M = m.m0;
else
    m1_mult();
    M = M + m.m0;
end

function m1_mult
    if size(m.m1, 3) == 1
        M = m.m1 * q;
    else
        M = zeros(size(m.m1, 1), size(m.m1, 3));
        for ii = 1 : size(m.m1, 3)
            M(:, ii) = m.m1(:, :, ii) * q;
        end
    end
end

end