function idx = detect_outliers(projections, theta, num_clusters,...
                               sigmaNoise, actual_classes, noisy_orientations,...
                               svector, max_angle_error, output_size)
    % % Agglomerative clustering.
    % Z = linkage(projections', 'centroid', 'euclidean');
    % num_clusters = inf;
    % for c=1.14:0.001:1.16
    %     idx = cluster(Z, 'cutoff', c);
    %     % idx = cluster(Z, 'Maxclust', 3);
        
    %     num_clus = size(unique(idx), 1);
        
    %     if (num_clus < num_clusters) && (num_clus > 300)
    %         num_clusters = num_clus;
    %         cutoff = c;
    %         curr_idx = idx;
    %     end
    % end
    
    % C = zeros(size(projections, 1), num_clusters);
    % idx = curr_idx;
    % for i=1:num_clusters
    %     projections_in_cluster = projections(:, idx == i);
    %     C(:, i) = mean(projections_in_cluster, 2);
    % end
    % clustered_projections = C';

    % K-means
    % [idx, C, ~, ~] = kmeans(projections', num_clusters,...
    %     'Replicates', 5, 'MaxIter', 1000);
    % clustered_projections = C';

    num_angles = size(projections, 2);
    epsilon = 30;
    projections = [projections flipud(projections)];

    W = pdist2(projections', projections', 'squaredeuclidean');
    W = exp(-W/(2*epsilon));

    % % Define the weight matrix.
    % W = zeros(num_angles, num_angles);
    % for i=1:num_angles
    %     for j=1:num_angles
    %         normie = norm(projections(:, i) - projections(:, j)).^2;
    %         W(i, j) = exp(-normie/(2*epsilon));
    %     end
    % end

    % Define diagonal matrix D.
    D = diag(sum(W, 2));
    W_tilda = D\W/D;
    D_tilda = diag(sum(W_tilda, 2));

    L = D_tilda - W_tilda;

    [V, ~] = eig(L, D_tilda);

    phi1 = -V(:, 2);
    phi2 = -V(:, 3);
    angles = (atan(phi1./phi2)*(180/pi) + 90)';
    angles = angles(1:num_angles);

    sorted_angles = sort(angles);

    idx = zeros(1, num_angles);
    for i=1:num_clusters
        if i == num_clusters
            div_value = inf;
        else
            div_value = sorted_angles(round((i+1)*num_angles/num_clusters -...
                (num_angles/num_clusters) + 1));
        end
        if i == 1
            prev_div_value = -inf;
        else
            prev_div_value = sorted_angles(round((i)*num_angles/num_clusters -...
                (num_angles/num_clusters) + 1));
        end
        idx((angles < div_value) & (angles >= prev_div_value)) = i;
    end

    C = zeros(size(projections, 1), num_clusters);
    idx = curr_idx;
    for i=1:num_clusters
        projections_in_cluster = projections(:, idx == i);
        C(:, i) = mean(projections_in_cluster, 2);
    end
    clustered_projections = C;

    % Calculate the angles of the clusters.
    clustered_angles = zeros(1, num_clusters);
    for i=1:num_clusters
        clustered_angles(i) = mean(theta(idx == i));
    end

    % Denoise the cluster centroids.
    sigmaNoise = sigmaNoise*num_clusters/size(theta, 2);
    measured_projections = denoise(clustered_projections, sigmaNoise, 145, 50);
    measured_projections = max(0, measured_projections);

    % If noisy_orientations is true...
    if noisy_orientations == 1
        initial_theta = clustered_angles +...
            randi([-max_angle_error max_angle_error], 1, num_clusters);
    end

    % Predict rough estimates of angle
    noisy_theta = ...
        ARPord(measured_projections, svector, sigmaNoise,...
        initial_theta, noisy_orientations, max_angle_error);

    no_shift = zeros(size(clustered_angles));
    reconstructed_image = ...
    reconstruct_image(clustered_projections, noisy_theta,...
        no_shift, output_size);
end