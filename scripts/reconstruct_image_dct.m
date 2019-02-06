function reconstructed_image = reconstruct_image_dct(...
    measured_projections, original_projections,...
    num_clusters, theta, sigmaNoise, num_theta, noisy_orientations,...
    max_shift_amplitude, svector, output_size, num_freq)
    
    % Specify the frequency domain approach.
    symmetry_method = 4;
    
    % Initial error between the projections.
    disp('**** L2-norm error between the original projections and measured projections ****');
    disp(norm(measured_projections - original_projections, 'fro'));
    disp('');

    disp('**** Cluster the projections ****');
    [clustered_projections, clustered_angles] =...
        cluster_projections(measured_projections, num_clusters, theta);

    % Save the original projections and mark the clustered projections as measured
    % projections.
    measured_projections = clustered_projections;
    theta = clustered_angles;

    % Calculate the new variance of noise in the projections.
    sigmaNoise = sigmaNoise*num_clusters/num_theta;

    % Update the number of clusters after filtering.
    num_clusters = size(clustered_projections, 2);

    % If orientations are noisy or completely unknown.
    if noisy_orientations == 1
        initial_theta = theta + randi([-1 1], 1, num_clusters);
    else
        initial_theta = randi([1 179], num_clusters, 1);
    end

    % The shifts estimated using the center of mass theorem.
    estimated_shifts = estimate_shifts(measured_projections, max_shift_amplitude); 
    original_shifts = estimated_shifts;

    % Error after shift correction.
    disp('**** L1-norm error between estimated shifts and actual shifts ****');
    disp(norm(estimated_shifts - original_shifts, 1));
    disp('');

    % Predict the angles and shifts using moment angle estimation.
    disp('**** Moment based estimation ****');
    [refined_projections, noisy_theta, projection_shifts] = ...
        SHARPord(measured_projections, svector, sigmaNoise, max_shift_amplitude,...
        estimated_shifts, initial_theta, noisy_orientations);
    disp('');
    projection_shifts = projection_shifts';
    noisy_theta = noisy_theta';

    disp('**** L1-norm error between the moment estimated shifts and the actual shifts ****');
    disp(norm(projection_shifts - original_shifts, 1));
    disp('');

    % Start iteration.
    projection_shifts = projection_shifts';
    noisy_theta = noisy_theta';

    better_theta = noisy_theta;
    better_shift = projection_shifts;

    % Define the gradient descent parameters.
    resolution_angle = 0.005;
    resolution_shift = 1;
    angle_amplitude = 3;
    shift_amplitude = 0;
    errors = [];
    alpha_rate =  0.001;
    beta_rate = 0.001;
    error_delta = 10;

    % Calculate the estimate of the image based on moment based solver.
    reconstructed_image = ...
        reconstruct_image(refined_projections, better_theta,...
            better_shift, output_size);

    % delta is the angle by which the image should be rotated
    % so that the axis of symmetry is horizontal.
    delta = estimate_axis_symmetry_alter(...
        reconstructed_image, output_size, 90, 0, symmetry_method);

    disp('**** Optimization error using symmetric prior ****');
    fprintf('\nIteration Error:            \n');
    for i=1:500

        % Use the projection prior
        gradient_vector = ...
            reconstructed_image -...
            reconstruct_image_dct_method(refined_projections, better_theta,...
            better_shift, output_size, num_freq);

        % Then use the symmetric prior.
        % For calculating the symmetry gradient we first rotate the image
        % to make the horizintal axis symmetrical.
        cropped_image = extract_circular_patch(reconstructed_image);

        rotated_image = imrotate(cropped_image, delta, 'bicubic', 'crop');
        symmetry_gradient_vector = symmetry_gradient_matrix(rotated_image);
        symmetry_gradient_vector = imrotate(symmetry_gradient_vector, -delta, 'bicubic', 'crop');

        % Now finally update the image.
        reconstructed_image = ...
            reconstructed_image - alpha_rate*gradient_vector - beta_rate*symmetry_gradient_vector; 

        % The optimization error for this iteration.
        function_error = calculate_optimization_error_symmetry_single_axis(refined_projections,...
            reconstructed_image, better_theta, better_shift, delta);

        % Display the error for ths iteration.
        fprintf('%6g \n', function_error); 

        % Do a brute force over all orientations and shifts.
        % Store the old estimated before calculating the new ones.
        [better_theta, better_shift] = best_orientation_estimate(...
            refined_projections, reconstructed_image, better_theta, better_shift,...
            angle_amplitude, shift_amplitude, resolution_angle,...
            resolution_shift);

        % Now update the symmetry axis.
        delta = estimate_axis_symmetry_alter(...
            reconstructed_image, output_size, error_delta, delta, symmetry_method);

        error_delta = max(1, error_delta - 0.5);
    end
end
