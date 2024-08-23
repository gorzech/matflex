function p = q_axis(alpha, axis)
%Q_AXIS Quaternions that rotate by alpha (in RAD) around axis (X-1, Y-2, Z-3)

e = zeros(3, 1);
e(axis) = sin(alpha / 2);
p = [cos(alpha / 2); e];