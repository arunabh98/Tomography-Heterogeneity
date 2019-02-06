image_size = 200;
P = imread('../../data/refs_017.png');
P = imresize(P, [image_size image_size]);
P = im2double(rgb2gray(P));

% Rotate the image.
P = imrotate(P, 20);

theta = 0:180;
shift = zeros(size(theta));
output_size = max(size(P));
projections = radon(P, theta);
theta = theta + randi([-5 5], 1, size(theta, 2));
delta_1_series = [];
delta_2_series = [];
delta_3_series = [];
delta_4_series = [];

%% Plot for noise
parfor i=1:30
    noise_fraction = 0.10 + i*0.01;
    noisy_projections = add_noise(projections, noise_fraction);
    delta_1 = estimate_axis_symmetry(noisy_projections, theta, shift,...
        output_size, 180, 1);
    delta_2 = estimate_axis_symmetry(noisy_projections, theta, shift,...
        output_size, 180, 2);
    delta_3 = estimate_axis_symmetry(noisy_projections, theta, shift,...
        output_size, 180, 3);
    delta_4 = estimate_axis_symmetry(noisy_projections, theta, shift,...
        output_size, 180, 4);
    delta_1_series = [delta_1_series delta_1];
    delta_2_series = [delta_2_series delta_2];
    delta_3_series = [delta_3_series delta_3];
    delta_4_series = [delta_4_series delta_4];
end

pdelta1 = process_deltas(delta_1_series);
pdelta2 = process_deltas(delta_2_series);
pdelta3 = process_deltas(delta_3_series);
pdelta4 = process_deltas(delta_4_series);
