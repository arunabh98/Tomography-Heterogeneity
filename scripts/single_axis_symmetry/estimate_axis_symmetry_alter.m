function delta = estimate_axis_symmetry_alter(...
	image_estimate, output_size, error_delta, prior_delta, method)

	% Extract the circular patch
	image_estimate = extract_circular_patch(image_estimate);
	
	min_diff = inf;
	for rot_ang=-error_delta:0.1:error_delta
		estimated_angle = mod(prior_delta + rot_ang, 180);
		rotated_image = imrotate(image_estimate, estimated_angle, 'bicubic', 'crop');

		diff_norm = amount_symmetric_horizontal(rotated_image, method);

		if diff_norm < min_diff
			min_diff = diff_norm;
			delta = estimated_angle;
		end
	end
end
