function delta = estimate_axis_symmetry(...
	projections, theta, shift, output_size, error_delta, method)

	min_diff = inf;
	for rot_ang=-error_delta:0.1:error_delta
		rotated_theta = mod(theta + rot_ang, 360);
		rotated_image = reconstruct_image(projections, rotated_theta,...
            shift, output_size);

		diff_norm = amount_symmetric_horizontal(rotated_image, method);

		if diff_norm < min_diff
			min_diff = diff_norm;
			delta = rot_ang;
		end
	end
end