function delta = estimate_double_axis_symmetry(...
	projections, theta, shift, output_size, error_delta, prior_delta, method)
	
	% Calculate the best axis and on the basis of that, the second axis.
	diff_rot_ang = zeros(1800, 1);
	c = 1;
	for rot_ang=0:0.1:179.9
		rotated_theta = mod(theta + rot_ang, 360);
		rotated_image = reconstruct_image(projections, rotated_theta,...
            shift, output_size);
		diff_rot_ang(c) = amount_symmetric_horizontal(rotated_image, method);
		c = c + 1;
	end

	% Sort the symmetry norms.
	[~, idx] = sort(diff_rot_ang);

	% Calculate the first estimate of our symmetry axis.
	first_estimates = zeros(2, 1);
	first_estimates(1) = (idx(1)/10) - 0.1;
	first_estimates(2) = mod(first_estimates(1) + 90, 180);

	% Start evaluation of the second estimates.
	second_estimates = zeros(2, 1);	
	% Find the next best estimate of symmetry - the second axis.
	for k=2:size(idx, 1)
		rot_angle_2 = (idx(k)/10) - 0.1;

		if abs(circ_dist(rot_angle_2, first_estimates(1))) > 50
			second_estimates(2) = rot_angle_2;
			second_estimates(1) = mod(rot_angle_2 - 90, 180);
            break;
		end
	end

	% Now we average the estimates.
	delta_1 = mod(first_estimates(1) + circ_dist(second_estimates(1), first_estimates(1))/2, 180);
	delta_2 = mod(first_estimates(2) + circ_dist(second_estimates(2), first_estimates(2))/2, 180);

	delta = delta_1;
	if prior_delta ~= -1
		if abs(circ_dist(delta_1, prior_delta)) < error_delta
			delta = delta_1;
		elseif abs(circ_dist(delta_2, prior_delta)) < error_delta
			delta = delta_2;
		elseif abs(circ_dist(delta_1, prior_delta)) < abs(circ_dist(delta_1, prior_delta))
			if circ_dist(delta_1, prior_delta) < 0
				delta = mod(prior_delta - error_delta, 180);
			elseif circ_dist(delta_1, prior_delta) > 0
				delta = mod(prior_delta + error_delta, 180);
			end			
		else
			if circ_dist(delta_2, prior_delta) < 0
				delta = mod(prior_delta - error_delta, 180);
			elseif circ_dist(delta_2, prior_delta) > 0
				delta = mod(prior_delta + error_delta, 180);
			end
		end
	end
end