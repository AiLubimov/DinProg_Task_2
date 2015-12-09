function S = gsqrtm(A)

[U, S, ~] = svd(A);
S = U * sqrt(S) * U';

end

