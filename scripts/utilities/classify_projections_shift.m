function [projection_classes, initial_incorrect, final_incorrect] = ...
    classify_projections_shift(projections, original_theta, original_class,...
    sigmaNoise, no_of_classes, filename)

    % Find the magnitude of fourier transform of each projection.
    f_projections = fft(projections);
    f_projections = fftshift(f_projections, 1);
    mag_projections = abs(f_projections);

    % Agglomerative clustering.
    Z = linkage(mag_projections', 'average', 'euclidean');
    num_clusters = inf;
    for c=1.14:0.001:1.16
        idx = cluster(Z, 'cutoff', c);
        
        num_clus = size(unique(idx), 1);
        
        if (num_clus < num_clusters) && (num_clus > 300)
            num_clusters = num_clus;
            curr_idx = idx;
        end
    end
    
    % Aggregate the clusters and statistics.
    C = zeros(size(projections, 1), num_clusters);
    class_clustered = zeros(1, num_clusters);
    theta_clustered = zeros(1, num_clusters);
    std_theta = zeros(1, num_clusters);
    std_class = zeros(1, num_clusters);
    zeroth_moment = zeros(1, num_clusters);
    idx = curr_idx;
    parfor i=1:num_clusters
        projections_in_cluster = projections(:, idx == i);
        C(:, i) = mean(projections_in_cluster, 2);
        class_clustered(i) = mode(original_class(idx == i));
        std_class(i) = std(original_class(idx == i));
        theta_clustered(i) = mean(original_theta(idx == i));
        std_theta(i) = std(original_theta(idx == i));
        zeroth_moment(i) = sum(C(:, i));
    end
    clustered_projections = C;

    % Denoise the projections.
    sigmaNoise = sigmaNoise*num_clusters/size(original_theta, 2);
    clustered_projections = denoise(clustered_projections, sigmaNoise, 145, 100);
    clustered_projections = max(0, clustered_projections);

    clustered_projections = [clustered_projections flipud(clustered_projections)];

    % Constants
    epsilon = 3e3;

    % Define the weight matrix.
    W = pdist2(clustered_projections', clustered_projections');
    W = W.^2;
    W = exp(-W/(2*epsilon));

    % Define diagonal matrix D.
    D = diag(sum(W, 2));
    W_tilda = D\W/D;
    D_tilda = diag(sum(W_tilda, 2));

    L = D_tilda - W_tilda;

    [V, ~] = eig(L, D_tilda);

    phi1 = -V(:, 2);
    phi2 = -V(:, 3);
    phi3 = -V(:, 4);
    coeff = [phi1 phi2 phi3];
    
    % Extract the three main eigenvectors.
    V = -V;
    V = V(1:num_clusters, :);
    coeff = coeff(1:num_clusters, :);
    phi1 = phi1(1:num_clusters);
    phi2 = phi2(1:num_clusters);
    phi3 = phi3(1:num_clusters);
    class_clustered = class_clustered';
    
    % Visualize the clusters.
    figure; scatter3(phi1, phi2, phi3, 10, class_clustered');
    savefig( strcat(filename, num2str(size(original_theta, 2)), '/original_clusters.fig'));

    % Initialize cluster classification using zeroth-order moment.
    if no_of_classes > 1
        % Initialize the cluster clases.
        [zeroth_moment_sorted, ~] = sort(zeroth_moment);
        thresholds = zeros(1, no_of_classes - 1);
        parfor i=1:(no_of_classes - 1)
            thresholds(i) = zeroth_moment_sorted(round(i*num_clusters/no_of_classes));
        end

        cluster_class = zeros(1, num_clusters);

        % The first class.
        cluster_class(zeroth_moment < thresholds(1)) =...
            mode(class_clustered(zeroth_moment < thresholds(1)));

        % Rest all classes.
        for i=2:no_of_classes-1
            cluster_class(zeroth_moment < thresholds(i) & zeroth_moment >=  thresholds(i-1)) =...
                mode(class_clustered(zeroth_moment < thresholds(i) & zeroth_moment >=  thresholds(i-1)));
        end

        % The last class.
        cluster_class(zeroth_moment >= thresholds(no_of_classes - 1)) =...
            mode(class_clustered(zeroth_moment >= thresholds(no_of_classes - 1)));
    else
        cluster_class = ones(1, new_num_clusters);
    end

    % Plot the initial classification of projections.
    figure; scatter3(phi1, phi2, phi3, 10, cluster_class);
    savefig(strcat(filename, num2str(size(original_theta, 2)), '/initial_classification.fig'));
    initial_incorrect = sum(class_clustered' ~= cluster_class);
    fprintf('Initial number of projections classified incorrectly: %d \r\n', initial_incorrect);

    % Nearest neighbour classification
    refined_cluster_class = zeros(1, num_clusters);
    neighbour_idx = knnsearch(coeff, coeff, 'K', 10, 'IncludeTies', true);
    for i = 1:size(coeff, 1)
        nearest_idx = cell2mat(neighbour_idx(i));
        nearest_classes = cluster_class(nearest_idx);
        refined_cluster_class(i) = mode(nearest_classes);
    end

    % Plot the refined cluster classes.
    figure; scatter3(phi1, phi2, phi3, 10, refined_cluster_class);
    savefig(strcat(filename, num2str(size(original_theta, 2)), '/final_classification.fig'));
    final_incorrect = sum(class_clustered' ~= refined_cluster_class);
    fprintf('Final number of projections classified incorrectly: %d \r\n', final_incorrect);
    cluster_class = refined_cluster_class;

    % Now assign the most likely class.
    for i=1:no_of_classes
        cluster_class(cluster_class == i) = no_of_classes + i;
    end
    for i=1:no_of_classes
        cluster_class(cluster_class == (no_of_classes + i)) = ...
            mode(class_clustered(cluster_class == (no_of_classes + i)));
    end

    projection_classes = zeros(1, size(projections, 2));
    % Now assign all projections their classes.
    for i=1:num_clusters
        class_cluster = cluster_class(i);
        projection_classes(idx == i) = class_cluster;
    end
end
