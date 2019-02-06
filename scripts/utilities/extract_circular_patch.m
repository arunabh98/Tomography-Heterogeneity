function cropped_image = extract_circular_patch(test_image)
	output_size = max(size(test_image));
	radius = output_size/2;
    [xx,yy] = ndgrid((1:output_size) - radius, (1:output_size) - radius);
    mask = (xx.^2 + yy.^2) < radius^2;
    cropped_image = test_image.*mask;
end
