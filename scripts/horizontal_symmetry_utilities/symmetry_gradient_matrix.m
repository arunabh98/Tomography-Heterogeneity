function symmetry_gradient = symmetry_gradient_matrix(test_matrix)
	half_size = floor(size(test_matrix, 1)/2);

	% Extract the first half.
	first_half = test_matrix(1:half_size, :);

	% Extract the second half. 
	flipped_matrix = flip(test_matrix, 1);
	second_half = flipped_matrix(1:half_size, :);

	% Calculate the difference between the first half and second half.
	diff_matrix = first_half - second_half;

	% Now create a negative symmetric test_matrix
	flipped_diff_matrix = -flip(diff_matrix, 1);

	if mod(size(test_matrix, 1), 2) == 1
		mid_row = zeros(1, size(test_matrix, 2));
		symmetry_gradient = [diff_matrix; mid_row; flipped_diff_matrix];
	else
		symmetry_gradient = [diff_matrix; flipped_diff_matrix];
	end

end