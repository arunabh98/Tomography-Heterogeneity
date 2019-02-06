function error_estimate = calculate_optimization_error_symmetry_double_axis(projections, image_estimate, theta, shift, delta)

	error_estimate = 0;
	parfor i=1:size(projections, 2)
		current_projection = projections(:, i);
		shifted_projection = circshift(current_projection, -shift(i));
		estimated_projection = radon(image_estimate, theta(i));

		% Add the error.
		error_estimate = error_estimate +...
			norm(shifted_projection - estimated_projection).^2;
	end

	% Now add the symmetry error.
	error_estimate = error_estimate +...
		amount_symmetric_optimization_function(image_estimate, delta);

	% Rotated image estimate.
	error_estimate = error_estimate +...
		amount_symmetric_optimization_function(image_estimate, delta + 90);
end