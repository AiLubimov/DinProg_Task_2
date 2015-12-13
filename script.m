% Initial data
clear
global t_0 t x_0 X_0 A_t B_t p_t P_t n_dim n_dir ;

k = 1;
m = 1;
b = 1;
x1 = 1;
x2 = 2;
x3 = 2;
v1 = -1;
v2 = 1;
v3 = 2;
U = 2;

t_0 = 0;
t = 8;
NSpan = 100;
tspan = linspace(t_0, t, NSpan);

 A_t = @(t) [0      0      0      1  0  0; 
             0      0      0      0  1  0; 
             0      0      0      0  0  1;
             -2*k/m k/m    k/m    -b 0  0;
             k/m    -2*k/m k/m    0  -b 0;
             k/m    k/m    -2*k/m 0  0  -b];
%A_t = @(t) eye(6);
%A_t = @(t) [5 0; 0 1];
 B_t = @(t) [0 0 0 0 0 0; 
             0 0 0 0 0 0;
             0 0 0 0 0 0;
             0 0 0 0 0 0;
             0 0 0 0 0 0;
             0 0 0 0 0 1];
%B_t = @(t) eye(2);        
p_t = @(t) [0 0 0 0 0 0]';
%p_t = @(t) [0 0]';
P_t = @(t) [0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 0;
            0 0 0 0 0 U];
%P_t = @(t) eye(2);
x_0 = [x1 x2 x3 v1 v2 v3]';
%x_0 = [0 0]';
X_0 = 0.001 * eye(6);

n_dir = 60;

n_dim = check_initial_data();

if n_dim == 0
    error('Wrong dimentions');
end

system = struct('A_t', A_t, 'B_t', B_t, 'P_t', P_t, 'p_t', p_t, ...
                'x_0', x_0, 'X_0', x_0, 'tspan', tspan, 'n_dim', n_dim);

% Gen directions
global L_0

basisMat = [1 0 0 0 0 0; 0 0 0 1 0 0]';

if rank(basisMat) ~= 2
    error('Linear dependence in basisMat');
end
if size(basisMat, 1) ~= n_dim
    error('Wrong dim of basisMat');
end

basisMat = normc(basisMat);
%basisMat(:, 1) = basisMat(:, 1) / norm(basisMat(:, 1));
%basisMat(:, 2) = basisMat(:, 2) / norm(basisMat(:, 2));

phi = linspace(0, 2 * pi, n_dir);
L_0 = [cos(phi); sin(phi)];
L_0 = basisMat * L_0;

%L_0 = 2 * (rand(n_dim, n_dir) - 0.5);
%L_0 = normc(L_0);

%%

clear
global t_0 t x_0 X_0 A_t B_t p_t P_t n_dim n_dir ;

t_0 = 0;
t = 1;
NSpan = 100;
tspan = linspace(t_0, t, NSpan);

A_t = @(t) [1 1 0;
            1 2 0;
            0 0 1];
B_t = @(t) eye(3);        
p_t = @(t) [0 0 0]';
P_t = @(t) eye(3);
x_0 = [0 0 0]';
X_0 = eye(3);

n_dir = 60;

n_dim = check_initial_data();

if n_dim == 0
    error('Wrong dimentions');
end

system = struct('A_t', A_t, 'B_t', B_t, 'P_t', P_t, 'p_t', p_t, ...
                'x_0', x_0, 'X_0', x_0, 'tspan', tspan, 'n_dim', n_dim);

% Gen directions
global L_0

basisMat = [1 0 0; 0 1 0]';

if rank(basisMat) ~= 2
    error('Linear dependence in basisMat');
end
if size(basisMat, 1) ~= n_dim
    error('Wrong dim of basisMat');
end

basisMat = normc(basisMat);
%basisMat(:, 1) = basisMat(:, 1) / norm(basisMat(:, 1));
%basisMat(:, 2) = basisMat(:, 2) / norm(basisMat(:, 2));

phi = linspace(0, 2 * pi, n_dir);
L_0 = [cos(phi); sin(phi)];
L_0 = basisMat * L_0;

%L_0 = 2 * (rand(n_dim, n_dir) - 0.5);
%L_0 = normc(L_0);

%%
% Calculate Tube
[T, EllCenCA, EllMatCA] = ReachTube(L_0, tspan);

%%
P = basisMat * inv(basisMat' * basisMat) * basisMat';
hold on
for i = 1 : n_dir
    q = P * EllCenCA{end, i};
    Q = P * EllMatCA{end, i} * P';
    X = getEllipsoidPoints(q, Q, 100, basisMat);
    %X = P * X;
    X = linsolve(basisMat, X);
    X = X';
    plot(X(:, 1), X(:, 2));
end
%%
% Draw tube. Static case.
figure
hold on
drawTube(T, EllCenCA, EllMatCA, basisMat, 'static');