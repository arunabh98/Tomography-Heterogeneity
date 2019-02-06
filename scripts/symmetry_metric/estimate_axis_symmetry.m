function delta = estimate_axis_symmetry(...
	projections, theta, shift, output_size, prior_delta, method)

	min_diff = inf;
	for rot_ang=-prior_delta:prior_delta
		rotated_theta = mod(theta + rot_ang, 360);
		rotated_image = iradon(projections, rotated_theta, output_size);

		diff_norm = amount_symmetric_horizontal(rotated_image, method);

		if diff_norm < min_diff
			min_diff = diff_norm;
			delta = rot_ang;
		end
	end
end