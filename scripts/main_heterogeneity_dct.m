close all;
clc;

% Increase the number of parpool workers.
% myCluster = parcluster('local');
% myCluster.NumWorkers = 28;
% saveProfile(myCluster); 
% parpool('local', 28);

addpath(genpath('../data'));
addpath(genpath('single_axis_symmetry'));
addpath(genpath('horizontal_symmetry_utilities'));
addpath(genpath('moment_based_estimation'));
addpath(genpath('noise_scripts'));
addpath(genpath('utilities'));

% Get the images of all the classes.
no_of_classes = 2;
image_size = 100;

P = zeros(image_size, image_size, no_of_classes);
parfor i=1:no_of_classes
    P(:, :, i) = read_process_image(strcat('proteins/protein_2/refs_00', num2str(i), '.png'), image_size);
end

% Constants.
non_uniform_distribution = 0;
sigmaNoiseFraction = 0.20;
if non_uniform_distribution == 0
    filename = ...
        strcat('../results/heterogeneity/num_class_', num2str(no_of_classes), '/', num2str(sigmaNoiseFraction*100), '_percent_noise/');
else
    filename = ...
        strcat('../results/heterogeneity/num_class_', num2str(no_of_classes), num2str(sigmaNoiseFraction*100), '_percent_noise/non_uniform_distribution/');
end
output_size = max(size(P(:, :, 1)));

% Experimentatal conditions.
symmetry_prior = 1;
noisy_orientations = 0;
symmetry_method = 4;
include_clustering = 1;
num_clusters = 540;
num_theta = 30000;
max_angle_error = 0;
max_shift_amplitude = 0;
num_freq = image_size/2;

% Create the folder to hold the results of the experiment.
mkdir(strcat(filename, num2str(num_theta), '/all_variables/'));
if include_clustering ~= 1
    theta_to_write = zeros(10, num_theta);
end

% The file which contains all the errors.
fileID = fopen(strcat(filename,...
    num2str(num_theta), '/result.txt'),'w');

% Write the original images.
parfor i=1:no_of_classes
    imwrite(P(:, :, i), strcat(filename, num2str(num_theta), '/original_image_', num2str(i), '.png'));
end

% Define ground truth angles and take the tomographic projection.
if non_uniform_distribution == 0
    theta = datasample(0:0.005:179, num_theta);
elseif non_uniform_distribution == 1
    theta = [datasample(0:0.005:20, num_theta/5)...
    datasample(40:0.005:60, num_theta/5)...
    datasample(80:0.005:120, 2*num_theta/5)...
    datasample(140:0.005:160, num_theta/5)];
end

% Generate the projections from all classes.
[projections, svector, original_class_of_projections] = ...
    get_projections(theta, P);

% [projections, svector] = radon(P,theta);
original_projections = projections;

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

disp('**** Initial classification of projections ****');
[projection_classes, ~, ~, ~] =...
    classify_projections(measured_projections, num_clusters, theta, original_class_of_projections,...
        sigmaNoise, no_of_classes);

disp('**** Initial classification performance ****')
fprintf('Number of projections classified correctly: %d',...
        sum(projection_classes ~= original_class_of_projections));

% No. of clusters to create while estimating the structure.
num_clusters = 180;

% Reconstruct the images.
for i=1:no_of_classes
    class_measured_projections = measured_projections(:, projection_classes == i);
    class_original_projections = original_projections(:, projection_classes == i);
    actual_classes = original_class_of_projections(projection_classes == i);
    class_theta = theta(projection_classes == i);
    class_num_theta = size(class_theta, 2);

    disp('**** Reconstructing image ****')
    reconstructed_image = reconstruct_image_dct(...
        class_measured_projections, class_original_projections,...
        num_clusters, class_theta, sigmaNoise, class_num_theta, noisy_orientations,...
        max_shift_amplitude, svector, output_size, num_freq);

    % Save the result.
    formatSpec = 'Final image rmse error: %4.2f \r\n';
    fprintf(fileID, formatSpec, calculate_rmse_error(reconstructed_image, P(:,:, i)));

    imwrite(reconstructed_image, ...
    strcat(filename, num2str(num_theta), '/reconstructed_image_',...
        num2str(i), '.png'));
end

% Save the important variables.
save(strcat(filename, num2str(num_theta), '/all_variables/all_variables.mat'),...
    '-regexp',...
    '^(?!(measured_projections|original_projections|projections)$).');

