function error_estimate  = calc_error(projections, kmax, svector, Ord, angles, class_estimate)

	% Seperate the projections and angles by their classes.
	[projections1, angles1] = extract_projections_and_angles_for_class(...
        projections, angles, class_estimate, 1);
    numkeep1 = size(projections1, 2);

    [projections2, angles2] = extract_projections_and_angles_for_class(...
        projections, angles, class_estimate, 2);
    numkeep2 = size(projections2, 2);

    [projections3, angles3] = extract_projections_and_angles_for_class(...
        projections, angles, class_estimate, 3);
    numkeep3 = size(projections3, 2);

    % Calculate the left hand side and the right hand side.
    PMord1 = ...
        assemblePMord(projections1, kmax, svector, Ord, numkeep1);
    PMord2 = ...
        assemblePMord(projections2, kmax, svector, Ord, numkeep2);
    PMord3 = ...
        assemblePMord(projections3, kmax, svector, Ord, numkeep3);

    A1 = assembleA(angles1, Ord);
    A2 = assembleA(angles2, Ord); 
    A3 = assembleA(angles3, Ord);

    IMestimated1 = A1 \ PMord1;
    IMestimated2 = A2 \ PMord2;
    IMestimated3 = A3 \ PMord3;

    % Calculate the error for each class.
    Ehlccvec1 = A1 * IMestimated1 - PMord1;
    Ehlccvec2 = A2 * IMestimated2 - PMord2; 
    Ehlccvec3 = A3 * IMestimated3 - PMord3;

    Ehlcc1 = norm(Ehlccvec1, 2);
    Ehlcc2 = norm(Ehlccvec2, 2);
    Ehlcc3 = norm(Ehlccvec3, 2);

    % Sum the error.
    error_estimate = Ehlcc1 + Ehlcc2 + Ehlcc3;
end