% Increase the number of parpool workers.
% parpool('local', 14)
warning('off', 'MATLAB:rankDeficientMatrix');

% Include the moment based estimation scripts and noise scripts.
addpath(genpath('../data'));
addpath(genpath('moment_based_estimation'));
addpath(genpath('noise_scripts'));
addpath(genpath('utilities'));

% Get the images
image_size = 100;
P1 = read_process_image('refs_008.png', image_size);
P2 = read_process_image('refs_009.png', image_size);
P3 = read_process_image('refs_010.png', image_size);

% Constants.
non_uniform_distribution = 0;
sigmaNoiseFraction = 0.40;
if non_uniform_distribution == 0
    filename = ...
        strcat('../results/heterogeneity/', num2str(sigmaNoiseFraction*100), '_percent_noise/');
else
    filename = ...
        strcat('../results/heterogeneity/', num2str(sigmaNoiseFraction*100), '_percent_noise/non_uniform_distribution/');
end
output_size = max(size(P1));

% Experimentatal conditions.
max_shift_amplitude = 0;
symmetry_prior = 1;
noisy_orientations = 1;
symmetry_method = 4;
include_clustering = 1;
num_clusters = 540;
num_theta = 30000;
max_angle_error = 0;

% Create the folder to hold the results of the experiment.
mkdir(strcat(filename, num2str(num_theta), '/all_variables/'));
if include_clustering ~= 1
    theta_to_write = zeros(10, num_theta);
end

% The file which contains all the errors.
fileID = fopen(strcat(filename,...
    num2str(num_theta), '/result.txt'),'w');

% Write the original images.
imwrite(P1, strcat(filename,...
    num2str(num_theta), '/original_image_1.png'));
imwrite(P2, strcat(filename,...
    num2str(num_theta), '/original_image_2.png'));
imwrite(P3, strcat(filename,...
    num2str(num_theta), '/original_image_3.png'));

% Define ground truth angles and take the tomographic projection.
if non_uniform_distribution == 0
    theta = datasample(0:0.005:179, num_theta);
elseif non_uniform_distribution == 1
    theta = [datasample(0:0.005:20, num_theta/5)...
    datasample(40:0.005:60, num_theta/5)...
    datasample(80:0.005:120, 2*num_theta/5)...
    datasample(140:0.005:160, num_theta/5)];
end

% Generate the projections including outliers of class 1 and 2.
[projections, svector, original_class_of_projections] = ...
    get_projections(theta, P1, P2, P3);

% [projections, svector] = radon(P,theta);
original_projections = projections;
original_shifts = zeros(size(theta));

% Shift each projection by an unknown amount.
parfor i=1:size(projections, 2)
    original_shifts(i) = ...
        randi([-max_shift_amplitude, max_shift_amplitude]);
    projections(:, i) = circshift(projections(:, i), original_shifts(i)); 
end

% Normalize s to a unit circle
smax = max(abs(svector));
svector = svector / smax;

% Define the length of the projection.
projection_length = size(projections, 1);

% Add noise to projections.
[measured_projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);

% Initial error between the projections.
disp('**** L2-norm error between the original projections and measured projections ****');
disp(norm(measured_projections - original_projections, 'fro'));
disp('');

if include_clustering == 1
    disp('**** Initial - Cluster the projections ****');
    [clustered_projections, clustered_angles, cluster_class, original_class_of_projections] =...
        cluster_projections(measured_projections, num_clusters, theta, original_class_of_projections,...
            sigmaNoise);

    % Save the original projections and mark the clustered projections as measured
    % projections.
    originally_measured_projections = measured_projections;
    measured_projections = clustered_projections;
    original_theta = theta;
    theta = clustered_angles;

    % Calculate the new variance of noise in the projections.
    sigmaNoise = sigmaNoise*num_clusters/num_theta;

    % Update the number of clusters after filtering.
    num_clusters = size(clustered_projections, 2);
    theta_to_write = zeros(10, num_clusters);
end

theta_to_write(1, :) = theta;
theta_to_write(2, :) = original_class_of_projections;

estimated_class_of_projections = cluster_class;
theta_to_write(3, :) = estimated_class_of_projections;

% % If orientations are noisy or completely unknown.
% if noisy_orientations == 1
%     if include_clustering == 1 
%         initial_theta = theta + randi([-max_angle_error max_angle_error], 1, num_clusters);
%     else 
%         initial_theta = theta + randi([-max_angle_error max_angle_error], 1, num_theta);
%     end
% else
%     if include_clustering == 1
%         initial_theta = randi([1 179], num_clusters, 1);
%     else 
%         initial_theta = randi([1 179], num_theta, 1);
%     end
% end

% % The shifts estimated using the center of mass theorem.
% estimated_shifts = estimate_shifts(measured_projections, max_shift_amplitude); 

% if include_clustering == 1
%     % Ignore shifts if we are clustering
%     original_shifts = estimated_shifts;
% end

% % Error after shift correction.
% disp('**** L1-norm error between estimated shifts and actual shifts ****');
% disp(norm(estimated_shifts - original_shifts, 1));
% disp('');     

% % Predict the angles and shifts using moment angle estimation.
% disp('**** Moment based estimation ****');
% measured_projections1 = measured_projections(:, estimated_class_of_projections == 1); 
% estimated_shifts1 = estimated_shifts(estimated_class_of_projections == 1);
% initial_theta1 = initial_theta(estimated_class_of_projections == 1);
% [refined_projections1, noisy_theta1, projection_shifts1] = ...
%     SHARPord(measured_projections1, svector, sigmaNoise, max_shift_amplitude,...
%     estimated_shifts1, initial_theta1, noisy_orientations, max_angle_error);

% measured_projections2 = measured_projections(:, estimated_class_of_projections == 2); 
% estimated_shifts2 = estimated_shifts(estimated_class_of_projections == 2);
% initial_theta2 = initial_theta(estimated_class_of_projections == 2);
% [refined_projections2, noisy_theta2, projection_shifts2] = ...
%     SHARPord(measured_projections2, svector, sigmaNoise, max_shift_amplitude,...
%     estimated_shifts2, initial_theta2, noisy_orientations, max_angle_error);

% measured_projections3 = measured_projections(:, estimated_class_of_projections == 3); 
% estimated_shifts3 = estimated_shifts(estimated_class_of_projections == 3);
% initial_theta3 = initial_theta(estimated_class_of_projections == 3);
% [refined_projections3, noisy_theta3, projection_shifts3] = ...
%     SHARPord(measured_projections3, svector, sigmaNoise, max_shift_amplitude,...
%     estimated_shifts3, initial_theta3, noisy_orientations, max_angle_error);
% disp('');

% better_shift1 = projection_shifts1';
% better_theta1 = noisy_theta1';

% better_shift2 = projection_shifts2';
% better_theta2 = noisy_theta2';

% better_shift3 = projection_shifts3';
% better_theta3 = noisy_theta3';

% better_theta = zeros(size(initial_theta));
% better_theta(estimated_class_of_projections == 1) = better_theta1;
% better_theta(estimated_class_of_projections == 2) = better_theta2;
% better_theta(estimated_class_of_projections == 3) = better_theta3;

% theta_to_write(4, :) = better_theta;

% % Calculate the estimate of the image based on moment based solver.
% reconstructed_image1 = ...
%     reconstruct_image(refined_projections1, better_theta1,...
%         better_shift1, output_size);

% reconstructed_image2 = ...
%     reconstruct_image(refined_projections2, better_theta2,...
%         better_shift2, output_size);

% reconstructed_image3 = ...
%     reconstruct_image(refined_projections3, better_theta3,...
%         better_shift3, output_size);

% % Define the gradient descent parameters.
% resolution_angle = 0.005;
% resolution_shift = 1;
% angle_amplitude = 3;
% shift_amplitude = 0;
% errors = [];
% alpha_rate = 0.001;
% beta_rate = 0.001;

% disp('**** Optimization error without using symmetric prior ****');
% fprintf('\nIteration Error:            \n');
% for i=1:100
%     prob_proj = calc_prob_for_proj(measured_projections, reconstructed_image1,...
%         reconstructed_image2, reconstructed_image3, better_theta);

%     % Use the projection prior
%     gradient_vector1 = ...
%         reconstructed_image1 -...
%         reconstruct_image_prob(measured_projections, better_theta, estimated_shifts, prob_proj(1, :), output_size);

%     % Now finally update the image.
%     reconstructed_image1 = ...
%         reconstructed_image1 - alpha_rate*gradient_vector1;

%     % Use the projection prior
%     gradient_vector2 = ...
%         reconstructed_image2 -...
%         reconstruct_image_prob(measured_projections, better_theta, estimated_shifts, prob_proj(2, :), output_size);

%     % Now finally update the image.
%     reconstructed_image2 = ...
%         reconstructed_image2 - alpha_rate*gradient_vector2;

%     % Use the projection prior
%     gradient_vector3 = ...
%         reconstructed_image3 -...
%         reconstruct_image_prob(measured_projections, better_theta, estimated_shifts, prob_proj(3, :), output_size);

%     % Now finally update the image.
%     reconstructed_image3 = ...
%         reconstructed_image3 - alpha_rate*gradient_vector3;

%      % Do a brute force over all orientations and shifts.
%     [better_theta1, better_shift1] = best_orientation_estimate(...
%         measured_projections, reconstructed_image1, better_theta, estimated_shifts,...
%         angle_amplitude, shift_amplitude, resolution_angle,...
%         resolution_shift);

%     [better_theta2, better_shift2] = best_orientation_estimate(...
%         measured_projections, reconstructed_image2, better_theta, estimated_shifts,...
%         angle_amplitude, shift_amplitude, resolution_angle,...
%         resolution_shift);

%     [better_theta3, better_shift3] = best_orientation_estimate(...
%         measured_projections, reconstructed_image3, better_theta, estimated_shifts,...
%         angle_amplitude, shift_amplitude, resolution_angle,...
%         resolution_shift);

%     better_theta = ...
%         prob_proj(1, :).*better_theta1 + prob_proj(2, :).*better_theta2 + prob_proj(3, :).*better_theta3;
%     better_shift = ...
%         prob_proj(1, :).*better_shift1 + prob_proj(2, :).*better_shift2 + prob_proj(3, :).*better_shift3;
% end

% [~, refined_class] = max(prob_proj);
% theta_to_write(5, :) = refined_class;
% theta_to_write(6, :) = better_theta;

% Write the thetas to csv file.
csvwrite(strcat(filename,...
    num2str(num_theta), '/thetas.csv'), theta_to_write);
fclose(fileID);

% % Save the important variables.
% save(strcat(filename, num2str(num_theta), '/all_variables/all_variables.mat'),...
%     '-regexp',...
%     '^(?!(clustered_projections|flipped_P1|flipped_P2|measured_projections|refined_projections|original_projections|originally_measured_projections|projections|first_quadrant|second_quadrant|third_quadrant|fourth_quadrant)$).');

