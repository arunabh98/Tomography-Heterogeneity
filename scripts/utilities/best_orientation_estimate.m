function [refined_theta, refined_shift] = best_orientation_estimate(...
	projections, image_estimate, theta_estimate, shift_estimate,...
	max_theta_amplitude, max_shift_amplitude, resolution_angle,...
	resolution_shift)
	
	% Choose the best theta and shift for each orientation.
	refined_theta = zeros(size(theta_estimate));
	refined_shift = zeros(size(shift_estimate));

	parfor i=1:size(projections, 2)

		% For each projection select the best orientation and shift.
		current_projection = projections(:, i);
		lowest_error = inf;

		for j=-max_shift_amplitude:resolution_shift:max_shift_amplitude

			% Iterate over all the possible shifts.
			current_shift = shift_estimate(i) + j;
			shifted_projection = circshift(current_projection, current_shift); 

			for k=-max_theta_amplitude:resolution_angle:max_theta_amplitude

				% Iterate over all possible orientations.
				current_theta = theta_estimate(i) + k;
				estimated_projection = radon(image_estimate, current_theta);

				if norm(estimated_projection - shifted_projection) < lowest_error
					lowest_error = norm(estimated_projection - shifted_projection);
					refined_theta(i) = current_theta;
					refined_shift(i) = -current_shift;
 				end
			end
		end
	end
end