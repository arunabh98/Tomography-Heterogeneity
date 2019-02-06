function err_t = err_for_all_angles(Pgiven, kmax, svector, Ord,...
    thetasestimated, max_limit, min_limit, i, PMord)

	range_size = round(max_limit - min_limit + 1);
	err_t = zeros(range_size, 1);

	parfor t = 1:range_size
        thetas_iter = thetasestimated;
        thetas_iter(i) = t + min_limit - 1;

        A = assembleA(thetas_iter, Ord);
        IMestimated = A \ PMord;
        E_tvec = A * IMestimated - PMord;
        err_t(t) = norm(E_tvec, 2);
    end
end