function s = sparsity(X)

X= X(:); % if matrix

%k = length(X);
%s = (sqrt(k) - (norm(X,1)/norm(X,2)) ) / ((sqrt(k) - 1));
s =  1 /  (norm(X,1));

end
