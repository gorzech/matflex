function plot_mode(q0, varargin)
%PLOT_MODE helper function to plot mode of the 3D beam with 6 DOFs
%   q0 - initial position of the beam
%   V - transformation matrix
%   idx - which mode to plot

if nargin < 3
    warning('Most likely there is nothing useful to be plotted')
end

q_value = 0.1; % What to set here?

plot_style = struct('n', 1, 'max', 4);
plot_style.styles = {'-o', '--d', '-.+', ':x'};

% Plot q0 as reference
q0_plot_style = '-*';
n_varargin = nargin - 1;
n_carg = 1;
if n_varargin >= n_carg
    if ischar(varargin{n_carg})
        q0_plot_style = varargin{n_carg};
        n_carg = n_carg + 1;        
    end    
end
if ~isempty(q0_plot_style)
    plot3(q0(1 : 6 : end), q0(2 : 6 : end), q0(3 : 6 : end), q0_plot_style, 'LineWidth', 1.5);
else
    clf
end

while n_varargin >= n_carg
    if isnumeric(varargin{n_carg})
        V = varargin{n_carg};
        n_carg = n_carg + 1;
    else
        error('Argument %d is expected to be numeric', n_carg)
    end
    if n_varargin < n_carg
        break;
    end
    if isnumeric(varargin{n_carg})
        idx = varargin{n_carg};
        if isscalar(idx)
            idx = [idx, idx];
        end
        n_carg = n_carg + 1;
    else
        error('Argument %d is expected to be scalar numeric', n_carg)
    end
    if n_varargin >= n_carg
        if isstring(varargin{n_carg})
            curr_plot_style = varargin{n_carg};
            n_carg = n_carg + 1;
        else
            curr_plot_style = get_next_style();
        end
    else
        curr_plot_style = get_next_style();
    end
    hold on
    if idx(1) > 0
        qm1 = q0 + V(:, idx(1)) * q_value;
        plot3(...
            qm1(1 : 6 : end), qm1(2 : 6 : end), qm1(3 : 6 : end), curr_plot_style, 'LineWidth', 1.5)
    end
    if idx(2) > 0
        qm2 = q0 - V(:, idx(2)) * q_value;
        plot3(...
            qm2(1 : 6 : end), qm2(2 : 6 : end), qm2(3 : 6 : end), curr_plot_style, 'LineWidth', 1.5)
    end    
    hold off
end

axis equal
view(30, 35)
grid on

    function s = get_next_style()
        s = plot_style.styles{plot_style.n};
        if plot_style.n >= 4
            plot_style.n = 1;
        else
            plot_style.n = plot_style.n + 1;
        end
    end
end