function error_mse = calculate_rmse_error(test_image, original_image)

    test_image = extract_circular_patch(test_image);
    original_image = extract_circular_patch(original_image);

    % Find the closest rotation to the image.
    dist_error = inf;
    flip_index = 0;
    for ang=-179:179
        rotated_test_image = imrotate(test_image, ang, 'bilinear', 'crop');
        flipped_image = flip(rotated_test_image, 1);
        if norm(rotated_test_image(:) - original_image(:), 2) < dist_error
            dist_error = norm(rotated_test_image(:) - original_image(:), 2);
            rot_ang = ang;
            flip_index = 0;
        end
        
        if norm(flipped_image(:) - original_image(:), 2) < dist_error
            dist_error = norm(flipped_image(:) - original_image(:), 2);
            rot_ang = ang;
            flip_index = 1;
        end
    end
    
    % Rotate the image.
    test_image = imrotate(test_image, rot_ang, 'crop');
    
    % If it is flipped.
    if flip_index == 1
        test_image = flip(test_image, 1);
    end

	% This function registers the image and then calculates the error.
	[optimizer,metric] = imregconfig('multimodal');
    optimizer.InitialRadius = optimizer.InitialRadius/3.5;
	test_image_registered = ...
		imregister(test_image, original_image, 'affine', optimizer, metric);

	error_mse = (norm(test_image_registered(:) - original_image(:))/norm(original_image(:)))*100.00;
end