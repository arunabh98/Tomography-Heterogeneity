
    % figure; gscatter(phi1, phi2, class_clustered);
    threshold = 0.02;
    max_size_cluster = 0;
    it = 0;
    while(it < 300)
%         random_index = find(class_clustered == 1);
        random_index = randi(num_clusters, round(num_clusters/(3)), 1);
        random_points = coeff(random_index, :);

        polynomial = polyfitn(random_points, ones(size(random_points, 1), 1), [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; 2 0 0; 0 2 0; 0 0 2]);        
        estimated_value = polyvaln(polynomial, coeff);

        size_cluster = sum(abs(estimated_value - 1) < threshold);

        if (size_cluster > max_size_cluster) && (size_cluster <= ceil(num_clusters/(2)))
            max_size_cluster = size_cluster;
            chosen_indices = coeff(abs(estimated_value - 1) < threshold, :);
            chosen_class = class_clustered(abs(estimated_value - 1) < threshold);
        end
        it = it + 1;
    end
    figure; scatter3(chosen_indices(:, 1), chosen_indices(:, 2), chosen_indices(:, 3), 10, chosen_class);