function n_dim = check_initial_data()
    global t_0 t x_0 X_0 A_t B_t p_t P_t n_dir;
    
    N = numel(x_0);
    
    n_dim = N;
    
    if t_0 >= t
        n_dim = 0;
        return;
    end
    
    if ~prod(size(A_t(t_0)) == [N N])
        n_dim = 0;
        return;
    end
    
    if ~prod(size(X_0) == [N N])
        n_dim = 0;
        return;
    end
    
    if numel(p_t(t_0)) ~= N
        n_dim = 0;
        return;
    end 
    
    if ~prod(size(P_t(t_0)) == [N N])
        n_dim = 0;
        return;
    end
    
    if ~prod(size(B_t(t_0)) == [N N])
        n_dim = 0;
        return;
    end
    
    if n_dir <= 0
        n_dim = 0;
        return;
    end
end

