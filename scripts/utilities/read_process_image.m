function processed_image = read_process_image(image_name, image_size)
	% Read the image.
	P = imread(strcat('data/', image_name));
	P = imresize(P, [image_size image_size]);
	P = im2double(rgb2gray(P));

	% Make sure the image is symmetric by taking an mirror image.
	% flipped_P2 = flip(P, 2);
	% P(:, (image_size/2) + 1:end) = flipped_P2(:, (image_size/2) + 1:end);

	% % Check if the image is symmetric.
	% left_half = P(:, 1:(image_size/2));
	% right_half = flip(P(:, (image_size/2) + 1:end), 2);
	% disp('**** Check if the image is symmetric ****');
	% disp(norm(left_half - right_half, 'fro'));
	% disp('');

	processed_image = P;

	% Rotate the image.
	% processed_image = imrotate(P, 90);
end