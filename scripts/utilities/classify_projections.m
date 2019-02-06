function [projection_classes, clustered_projections, clustered_angles, cluster_class] = ...
    classify_projections(projections, num_clusters, original_theta, original_class,...
    sigmaNoise, no_of_classes)

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

    if no_of_classes > 1
        % Initialize the cluster clases.
        [zeroth_moment_sorted, ~] = sort(zeroth_moment);
        thresholds = zeros(1, no_of_classes - 1);
        parfor i=1:(no_of_classes - 1)
            thresholds(i) = zeroth_moment_sorted(round(i*new_num_clusters/no_of_classes));
        end

        cluster_class = zeros(1, new_num_clusters);

        % The first class.
        cluster_class(zeroth_moment < thresholds(1)) =...
            mode(original_cluster_class(zeroth_moment < thresholds(1)));

        % Rest all classes.
        for i=2:no_of_classes-1
            cluster_class(zeroth_moment < thresholds(i) & zeroth_moment >=  thresholds(i-1)) =...
                mode(original_cluster_class(zeroth_moment < thresholds(i) & zeroth_moment >=  thresholds(i-1)));
        end

        % The last class.
        cluster_class(zeroth_moment >= thresholds(no_of_classes - 1)) =...
            mode(original_cluster_class(zeroth_moment >= thresholds(no_of_classes - 1)));
    else
        cluster_class = ones(1, new_num_clusters);
    end;

    % Analyze cluster purity.
    original_cluster_purity = zeros(1, num_clusters);

    parfor ang = 1:num_clusters
        classProjs = original_class(idx == ang);

        frequent_original_class = mode(classProjs);
        original_cluster_purity(ang) =...
            (sum(classProjs == frequent_original_class)/size(classProjs, 2))*100;
    end
    
    disp('');
    disp(mean(original_cluster_purity, 2));

    projection_classes = zeros(1, size(projections, 2));
    % Now assign all projections their classes.
    for i=1:num_clusters
        class_cluster = cluster_class(i);
        projection_classes(idx == i) = class_cluster;
    end
end