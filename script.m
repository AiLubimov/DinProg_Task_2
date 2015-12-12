% Initial data
clear
global t_0 t x_0 X_0 A_t B_t p_t P_t n_dim n_dir ;

t_0 = 0;
t = 1;
NSpan = 50;
tspan = linspace(t_0, t, NSpan);

A_t = @(t) [1, 0, 0;
            0, 1, 0;
            0, 0, 1];
B_t = @(t) [1, 0, 0;
            0, 1, 0;
            0  0, 1];
        
p_t = @(t) [0 0 0]';
P_t = @(t) [1 0 0;
            0 1 0;
            0 0 1];

x_0 = [0 0 0]';
X_0 = eye(3);

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
L_0 = [cos(phi); sin(phi)];

basisMat = [1 4 1; 2 1 6]';

if rank(basisMat) ~= 2
    error('Linear dependence in basisMat');
end
if size(basisMat, 1) ~= n_dim
    error('Wrong dim of basisMat');
end

basisMat = normc(basisMat);
%basisMat(:, 1) = basisMat(:, 1) / norm(basisMat(:, 1));
%basisMat(:, 2) = basisMat(:, 2) / norm(basisMat(:, 2));

L_0 = basisMat * L_0;

%%
% Calculate Tube
[T, EllCenCA, EllMatCA] = ReachTube(L_0, tspan);

%%
% Draw tube. Static case.
figure
hold on
drawTube(T, EllCenCA, EllMatCA, basisMat, 'dynamic', A_t);