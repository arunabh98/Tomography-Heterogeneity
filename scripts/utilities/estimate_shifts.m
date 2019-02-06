function estimated_shifts = estimate_shifts(projections, max_shift_amplitude)

    projection_length = size(projections, 1);
    estimated_shifts = zeros(1, size(projections, 2));

    % Estimate the shifts by keeping the center of mass in the center.
    parfor i=1:size(projections, 2)
        current_projection = projections(:, i);

        % Calculate center of mass of the projection.
        tot_mass = sum(current_projection(:));
        [ii, ~] = ...
            ndgrid(1:size(current_projection,1),1:size(current_projection,2));
        center_of_mass = sum(ii(:).*current_projection(:))/tot_mass;

        % Calculate the shift amount and limit it by what we
        % know to be the maximum.
        shift_amount = round(((projection_length + 1)/2) - center_of_mass);
        if shift_amount > max_shift_amplitude
            shift_amount  = max_shift_amplitude;
        elseif shift_amount < -max_shift_amplitude
            shift_amount  = -max_shift_amplitude;
        end
        estimated_shifts(i) = -shift_amount;
    end
end