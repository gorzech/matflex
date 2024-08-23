function animate_mode(q0, V, idx)
%ANIMATE_MODE helper function to animate modes of the 3D beam with 6 DOFs
%   q0 - initial position of the beam
%   V - transformation matrix
%   idx - which mode to animate

q_value = 0.1; % What to set here?
qm1 = q0 - V(:, idx) * q_value;
qm2 = q0 + V(:, idx) * q_value;

q_all = [q0, qm1, qm2];
plot_xlim = [min(min(q_all(1:6:end, :))), max(max(q_all(1:6:end, :)))];
plot_ylim = [min(min(q_all(2:6:end, :))), max(max(q_all(2:6:end, :)))];
plot_zlim = [min(min(q_all(3:6:end, :))), max(max(q_all(3:6:end, :)))];

% plot in 3d q0 and qm
p = plot3(q0(1 : 6 : end), q0(2 : 6 : end), q0(3 : 6 : end), '-*', ...
    qm1(1 : 6 : end), qm1(2 : 6 : end), qm1(3 : 6 : end), '-o');

axis equal
view(30, 35)
grid on
xlim(plot_xlim)
ylim(plot_ylim + [-0.1, 0.1])
zlim(plot_zlim + [-0.1, 0.1])

% Animate
dq_values = linspace(-q_value, q_value, 10);
dq_values = [dq_values, dq_values(end - 1 : -1 : 1)];
for dq = dq_values
    qm = q0 + V(:, idx) * dq;
    set(p(2), 'XData', qm(1 : 6 : end));
    set(p(2), 'YData', qm(2 : 6 : end));
    set(p(2), 'ZData', qm(3 : 6 : end));
    pause(0.1);
end

