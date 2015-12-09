% Initial data
global t_0 t x_0 X_0 A_t B_t p_t P_t n_dim n_dir P;

t_0 = 0;
t = 1.5;
NSpan = 100;
tspan = linspace(t_0, t, NSpan);

A_t = @(t) [0, 0;
            0, 1-t];
B_t = @(t) [0, 0;
            t, 0];
        
p_t = @(t) [0 0]';
P_t = @(t) [1 1;
            2 2];

x_0 = [0 0]';
X_0 = eye(2);

n_dir = 30;

n_dim = check_initial_data();

if n_dim == 0
    error('Wrong dimentions');
end

system = struct('A_t', A_t, 'B_t', B_t, 'P_t', P_t, 'p_t', p_t, ...
                'x_0', x_0, 'X_0', x_0, 'tspan', tspan, 'n_dim', n_dim);

% Gen directions
global L_0


phi = linspace(0, 2 * pi, n_dir);
L_0 = [sin(phi); cos(phi)];

basisMat = [1 0; 0 1]';
L_0 = basisMat * L_0;

%L_0 = (rand(n_dim, n_dir) - 0.5) * 2;
%L_0(:, i) = L_0(:, i) / norm(L_0(:, i));
%L_0(3 :4, :) = zeros(2, n_dir);
%P = basisMat * inv(basisMat' * basisMat) * basisMat';
%%
% Calculate Tube
[T, EllCenCA, EllMatCA] = ReachTube(L_0, tspan);

%%
% Draw tube. Static case.

drawTube_static(T, EllCenCA, EllMatCA, basisMat);