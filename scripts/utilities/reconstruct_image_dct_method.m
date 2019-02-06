function image_estimate_dct = reconstruct_image_dct_method(...
	projections, theta_estimate, shift_estimate, output_size, num_freq)

	shifted_projections = correct_projection_shifts(projections, shift_estimate);
	image_estimate = iradon(shifted_projections, theta_estimate, output_size, 'Cosine');
	D = dctmtx(output_size);

    coeff = D'*image_estimate*D;

    image_estimate_dct = zeros(output_size, output_size);
    for i = 1:num_freq
        for j = 1:num_freq
            basis_matrix = D(i,:)' * D(j,:);
            image_estimate_dct = image_estimate_dct + coeff(i, j)*basis_matrix;
        end
    end
end
