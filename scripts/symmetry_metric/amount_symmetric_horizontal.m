function symmetry_norm = amount_symmetric_horizontal(test_image, method)

	if method == 1
		%%%%% The L2 norm estimate. %%%%%
		% Extract the first half.
		height = size(test_image, 2);
		first_half_length = floor(height/2);
		first_half = test_image(1:first_half_length, :);
		
		% Extract the second half.
		flipped_image = flip(test_image, 1);
		second_half = flipped_image(1:first_half_length, :);

		% Return the difference.
		symmetry_norm = norm(second_half - first_half).^2;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif method == 2
		%%%% The L1 norm estimate. %%%%%
		% Extract the first half.
		height = size(test_image, 2);
		first_half_length = floor(height/2);
		first_half = test_image(1:first_half_length, :);
		
		% Extract the second half.
		flipped_image = flip(test_image, 1);
		second_half = flipped_image(1:first_half_length, :);

		% Return the difference.
		symmetry_norm = norm(second_half - first_half, 1);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif method == 3
		%%%%% The L0.5 norm estimate. %%%%%
		% Extract the first half.
		height = size(test_image, 2);
		first_half_length = floor(height/2);
		first_half = test_image(1:first_half_length, :);
		
		% Extract the second half.
		flipped_image = flip(test_image, 1);
		second_half = flipped_image(1:first_half_length, :);

		% Return the difference.
        diff_mat = (second_half - first_half).^(0.5);
        diff_sum = sum(diff_mat(:));
		symmetry_norm = diff_sum^2;
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	elseif method == 4
		symmetry_norm = measure_symmetry(test_image);
	end


end