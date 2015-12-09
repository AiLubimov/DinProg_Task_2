function drawTube_static(T, EllCenCA, EllMatCA, basisMat)

    P = basisMat * inv(basisMat' * basisMat) * basisMat';
    ConvHullCA = cell(numel(T), 1);
    %DirCA = cell(numel(T), 1);
    %DirArr = zeros(n_dir, 2);
    n_dir = size(EllCenCA, 2);
    
    for k = 1 : numel(T)
        X = [];
        for i = 1 : n_dir
            EllCen = P * EllCenCA{k, i};
            EllMat = P * EllMatCA{k, i} * P';
        
            X = [X getEllipsoidPoints(EllCen, EllMat, 100)];
        
            ind = convhull(X(1, :)', X(2, :)');
            X = [X(1, ind); X(2, ind)];
            %x = EllCen + EllMat * L_0(:, i) / dot(L_0(:, i), EllMat * L_0(:, i)) ^ 0.5;
            %DirArr(i, :) = x;
        end
        ConvHullCA{k} = X;
        %DirCA{k} = DirArr;
    end

    max = 0;

    for i = 1 : numel(ConvHullCA)
        if max < size(ConvHullCA{i}, 2)
            max = size(ConvHullCA{i}, 2);
        end
    end

    for i = 1 : numel(ConvHullCA)
        convHull = ConvHullCA{i};
        if size(convHull, 2) < max
            convHull = [interp1(linspace(0, 1, size(convHull,2)), convHull(1, :), linspace(0, 1, max));
            interp1(linspace(0, 1, size(convHull,2)), convHull(2, :), linspace(0, 1, max))];
        end
        ConvHullCA{i} = convHull;
    end

    Z = zeros(max, numel(T));
    X = zeros(max, numel(T));
    Y = zeros(max, numel(T));

    for i = 1 : numel(T) 
        convHull = ConvHullCA{i};
        X(:, i) = convHull(1, :);
        Y(:, i) = convHull(2, :);
        Z(:, i) = T(i) * ones(max, 1);
    end

    surf(X, Y, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

end

