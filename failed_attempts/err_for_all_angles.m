function err_t = err_for_all_angles(shiftedPgiven, kmax, svector, Ord,...
    thetasestimated, class_iter, max_limit, min_limit, i)

	range_size = round(max_limit - min_limit + 1);
	err_t = zeros(range_size, 1);

	parfor t = 1:range_size
        warning('off', 'MATLAB:rankDeficientMatrix');
        thetas_iter = thetasestimated;
        thetas_iter(i) = t + min_limit - 1;

        err_t(t)  = calc_error(...
            shiftedPgiven, kmax, svector, Ord,...
            thetas_iter, class_iter);
    end
end