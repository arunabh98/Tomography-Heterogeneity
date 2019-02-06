% DBSCAN Clustering algorithm.
epsilon=0.07;
MinPts=10;
cluster_class=DBSCAN(coeff,epsilon,MinPts);

% Kernel K-means clustering.
init = cluster_class;
cluster_class = knKmeans(coeff', init, @knPoly);

% Visualize the clusters.
no_of_classes = size(unique(cluster_class), 1);
figure; scatter3(phi1, phi2, phi3, 10, cluster_class);

% it = 0;
    % % best_clustering
    % deviation = inf;
    % while(it < 500)
    %     random_index = randperm(num_clusters, round(num_clusters/(rand() + 1)));
    %     random_points = coeff(random_index, :);

    %     % Cluster the random points.
    %     Z = linkage(random_points, 'single', 'euclidean');
    %     c = cluster(Z, 'Maxclust', no_of_classes);
        
    %     min_num_cluster = min(hist(c, unique(c)));

    %     if abs(min_num_cluster - (num_clusters/no_of_classes)) < deviation
    %         deviation = abs(min_num_cluster - (num_clusters/no_of_classes));
    %         cluster_class = zeros(1, num_clusters);
    %         cluster_class(random_index) = c;
    %         best_idx = c;
    %         best_random_points = random_points;
    %         best_random_index = random_index;
    %     end
        
    %     it = it + 1;
    % end
    
    % if deviation > 150
    %     disp('Clusters cannot be automatically identified.');
    %     return;
    % end

    % % Points not identified.
    % index = find(cluster_class == 0);
    % remaining_points = coeff(index, :);
    % % Identify the remaining points.
    % [~, I] = pdist2(best_random_points, remaining_points, 'euclidean', 'Smallest', 1);
    
    % cluster_class(index) = best_idx(I);