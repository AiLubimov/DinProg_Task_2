function [T, EllCenCA, EllMatCA] = ReachTube(L_0, tspan)

global x_0 X_0 n_dim n_dir
%global dotl0X0l02 l_0

EllCenCA = cell(numel(tspan), n_dir);
EllMatCA = cell(numel(tspan), n_dir);

X_0_Vec = reshape(gsqrtm(X_0), [n_dim * n_dim 1]);

for i = 1 : n_dir
    
    l_0 = L_0(:, i);
    dotl0X0l02 = dot(l_0, X_0 * l_0) ^ 0.5;
    
    [T, Ell] = ode45(@odefun, tspan, [l_0; x_0; X_0_Vec]);
    
    for k = 1 : numel(T)
        EllCenCA{k, i} = Ell(k, n_dim + 1 : 2 * n_dim)';
        EllMatS = reshape(Ell(k, 2 * n_dim + 1 : end)',[n_dim n_dim]);
        EllMatCA{k, i} = EllMatS' * EllMatS;
    end
end
    
function dy = odefun(t, y)

    global A_t B_t p_t P_t

    At = A_t(t); Bt = B_t(t); pt = p_t(t); Pt = P_t(t);

    P_b = Bt * Pt * Bt'; 
    P_b_s = gsqrtm(P_b);

    l = y(1 : n_dim);
    q_ = y(n_dim + 1 : 2 * n_dim);
    Qs = reshape(y(2 * n_dim + 1 : end), [n_dim n_dim]);

    %-Calc S------------------------------
    %global dotl0X0l02 l_0;

    if abs(dotl0X0l02) > eps
        lamt = dot(l, P_b * l)^0.5 / dotl0X0l02;
    else 
        lamt = 0;
    end

    E = eye(n_dim);

    if abs(lamt) > eps
    
        v_2 = P_b_s * l;
        v_1 = lamt * sqrtm(X_0) * l_0; 

        nv_1 = v_1 / norm(v_1);
        nv_2 = v_2 / norm(v_2);
        c = dot(nv_1, nv_2);
        s = sqrt(1 - (c^2 - 10^(-15)));
        q_1 = nv_1;
        
        if abs(s) > eps
            q_2 = (nv_2 - c * nv_1) / s;
        else
            q_2 = zeros(n_dim, 1);
        end

        Q_1 = [q_1 q_2];

        S = E + Q_1 * ([c s; -s c] - [1 0; 0 1]) * Q_1';
    else
        S = E;
    end
    %--------------------------------------
    
    dl = -At' * l;
    dl = dl / norm(dl);
    dq_ = At * q_ + Bt * pt;
    dQs = S * P_b_s + Qs * At';

    dy = [dl; dq_; reshape(dQs, [n_dim * n_dim, 1])];
end

end

