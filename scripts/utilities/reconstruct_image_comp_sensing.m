function image_estimate = reconstruct_image_comp_sensing(...
	projections, theta_estimate, shift_estimate, output_size, height, width, projection_length)
	
	% Define the basis matrix.
	D = dctmtx(height);

	% Correct the shifts.
	shifted_projections = correct_projection_shifts(projections, shift_estimate);

	% Reconstruct the image using compressed sensing. 

	% Define the parameters.
	n = height*width;
    m = projection_length*size(theta_estimate, 2);
    lambda  = 0.1;
	rel_tol = 200;

    % Set up the variables.
	y = shifted_projections(:);
	A = radonTransform(...
        	theta_estimate', width, height, output_size, projection_length);
    At = A';

	%run the l1-regularized least squares solver
    [reconstructed_image, ~]= ...
        l1_ls(A,At,m,n,y,lambda,rel_tol,true);

    % Reconstruct the image.
    image_estimate = reshape(reconstructed_image, [height, width]);
    image_estimate = D'*image_estimate;
    image_estimate(image_estimate < 0) = 0;
end