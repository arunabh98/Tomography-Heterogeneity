function symmetry_norm = amount_symmetric_optimization_function(test_image, delta)

	% Extract the patch which matters.
	test_image = imrotate(test_image, delta, 'bicubic', 'crop');
	test_image = extract_circular_patch(test_image);

	%%%%% The L2 norm estimate. %%%%%
	% Extract the first half.
	height = size(test_image, 2);
	first_half_length = floor(height/2);
	first_half = test_image(1:first_half_length, :);
	
	% Extract the second half.
	flipped_image = flip(test_image, 1);
	second_half = flipped_image(1:first_half_length, :);

	% Return the difference.
	symmetry_norm = norm(second_half(:) - first_half(:)).^2;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end