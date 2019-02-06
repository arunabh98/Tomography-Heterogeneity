function image_estimate = reconstruct_image_prob(...
	projections, theta_estimate, shift_estimate, prob_estimate, output_size)

	shifted_projections = correct_projection_shifts(projections, shift_estimate);
	prob_proj = repmat(prob_estimate, size(projections, 1), 1);
	prob_weighted_proj = shifted_projections.*prob_proj;
	image_estimate = iradon(prob_weighted_proj, theta_estimate, output_size, 'Cosine');
end