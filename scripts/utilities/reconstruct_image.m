function image_estimate = reconstruct_image(...
	projections, theta_estimate, shift_estimate, output_size)

	shifted_projections = correct_projection_shifts(projections, shift_estimate);
	image_estimate = iradon(shifted_projections, theta_estimate, output_size, 'Cosine');
end