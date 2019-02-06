function shifted_projections = correct_projection_shifts(...
	projections, estimated_shifts)

	shifted_projections = zeros(size(projections));

	parfor i=1:size(projections, 2)
		shifted_projections(:, i) = ...
			circshift(projections(:, i), -estimated_shifts(i)); 
	end

end