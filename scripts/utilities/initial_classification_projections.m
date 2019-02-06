function [projection_classes, clustered_projections, clustered_angles, cluster_class] =...
    initial_classification_projections(projections, num_clusters,...
    original_theta, original_class, sigmaNoise)
    % K-means
    [idx, clustered_projections, ~, ~] = kmeans(projections', num_clusters,...
        'distance', 'sqeuclidean', 'Replicates', 5, 'MaxIter', 1000);
    clustered_projections = clustered_projections';

    % Denoise the clusters.
    % Calculate the new variance of noise in the projections.
    sigmaNoise = sigmaNoise*num_clusters/size(original_theta, 2);

    clustered_projections = denoise(clustered_projections, sigmaNoise, 145, 50);
    clustered_projections = max(0, clustered_projections);

    % Calculate the angles and class of projections in each cluster.
    new_num_clusters = size(clustered_projections, 2);
    clustered_angles = zeros(1, new_num_clusters);
    original_cluster_class = zeros(1, new_num_clusters);
    zeroth_moment = zeros(1, new_num_clusters);
    parfor i=1:new_num_clusters
        clustered_angles(i) = mean(original_theta(idx == i));
        original_cluster_class(i) = mode(original_class(idx == i));
        zeroth_moment(i) = sum(clustered_projections(:, i));
    end

    % Initialize the cluster clases.
    [zeroth_moment_sorted, ~] = sort(zeroth_moment);
    first_threshold = zeroth_moment_sorted(round(new_num_clusters/3));
    second_threshold = zeroth_moment_sorted(round(2*new_num_clusters/3));
    cluster_class = zeros(1, new_num_clusters);
    cluster_class(zeroth_moment < first_threshold) =...
        mode(original_cluster_class(zeroth_moment < first_threshold));
    cluster_class(zeroth_moment < second_threshold & zeroth_moment >= first_threshold) =...
        mode(original_cluster_class(zeroth_moment < second_threshold & zeroth_moment >= first_threshold));
    cluster_class(zeroth_moment >= second_threshold) =...
        mode(original_cluster_class(zeroth_moment >= second_threshold));

    % Analyze cluster purity.
    original_cluster_purity = zeros(1, num_clusters);

    parfor ang = 1:num_clusters
        classProjs = original_class(idx == ang);

        frequent_original_class = mode(classProjs);
        original_cluster_purity(ang) =...
            (sum(classProjs == frequent_original_class)/size(classProjs, 2))*100;
    end
    
    disp(mean(original_cluster_purity, 2));

    projection_classes = zeros(1, size(projections, 2));
    % Now assign all projections their classes.
    for i=1:num_clusters
        class_cluster = cluster_class(i);
        projection_classes(idx == i) = class_cluster;
    end

    % Display number of wrong predictions.
    disp(sum(projection_classes ~= original_class));
end