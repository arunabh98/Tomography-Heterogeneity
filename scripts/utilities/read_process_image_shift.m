function processed_image = read_process_image_shift(image_name, image_size)
	% Read the image.
	P = imread(strcat('data/', image_name));
	P = imresize(P, [image_size image_size]);
	P = im2double(rgb2gray(P));

	% Pad the image with a fixed boundary of 5 pixels.
	P = padarray(P, [3, 3], 0.0);
	processed_image = P;

	% Rotate the image.
	% processed_image = imrotate(P, 90);
end