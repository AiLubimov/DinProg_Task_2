function X = getEllipsoidPoints(q, Q, N, basisMat) 
    
    phiphi = linspace(0, 2 * pi, N);

    x(1, :) = cos(phiphi);
    x(2, :) = sin(phiphi);
    
    if ~isempty(basisMat)
        x = basisMat * x;
    end
    
    [U, S, V] = svd(Q);
    X = repmat(q, 1, N) + U * sqrt(S) * V' * x;
    %plot(X(1, :), X(2, :), 'r')
end