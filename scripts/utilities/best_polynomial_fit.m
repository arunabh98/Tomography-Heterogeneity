function best_polynomial = best_polynomial_fit(points)
	most_points = 0;
	for it=1:500
		random_index = randperm(size(points, 1), round((3/4)*size(points, 1)));
		random_points = points(random_index, :);

		% polynomial = ...
        %     fit(random_points(:, 1:2), random_points(:, 3), 'loess',...
        %     'Robust', 'Bisquare', 'Span', 0.4);
        %  estimated_value = feval(polynomial, points(:, 1:2));

        polynomial = polyfitn(random_points(:, 1:2), random_points(:, 3), 2);
		estimated_value = polyvaln(polynomial, points(:, 1:2));

		diff_estimated = abs(estimated_value - points(:, 3));
		no_points_inside = sum(diff_estimated < 0.03);

		if no_points_inside > most_points
			best_polynomial = polynomial;
			most_points = no_points_inside;
		end
	end
end