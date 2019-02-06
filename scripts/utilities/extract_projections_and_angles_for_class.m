function [projections, theta] = extract_projections_and_angles_for_class(...
	all_projections, all_angles, class_projections, class_required)
		
	projections = all_projections(:, find(class_projections == class_required));
	theta = all_angles(:, find(class_projections == class_required));
end