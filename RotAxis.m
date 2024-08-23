function res = RotAxis(alpha, axis)
%ROTAXIS rotate by alpha (in RAD) around axis (X-1, Y-2, Z-3)
    
    ca = cos(alpha);
    sa = sin(alpha);
    res = zeros(3)*ca;

    if axis == 1
        res(2, 2) = ca;
        res(3, 3) = ca;
        res(2, 3) = -sa;
        res(3, 2) = sa;
    elseif axis == 2
        res(1, 1) = ca;
        res(3, 3) = ca;
        res(1, 3) = sa;
        res(3, 1) = -sa;
    elseif axis == 3
        res(1, 1) = ca;
        res(1, 2) = -sa;
        res(2, 1) = sa;
        res(2, 2) = ca;
    else
        error('Wrong axis specification in call to RotAxis')
    end
    res(axis, axis) = 1;
end