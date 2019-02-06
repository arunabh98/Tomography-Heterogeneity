function [clustered_projections, clustered_angles] = ...
    cluster_projections(projections, num_clusters, original_theta, outlier_mode, outlier_percentage, outlier_indices)
    % Cluster the projections.
    [idx, C, ~, ~] = kmeans(projections', num_clusters,...
        'distance', 'cityblock');
    if outlier_mode ~= 1
	    clustered_angles = zeros(1, num_clusters);
	    for i=1:num_clusters
	        clustered_angles(i) = mean(original_theta(find(idx == i)));
	    end
	    clustered_projections = C';
	else
		% Make the Distance Matrix
		D = zeros(size(projections, 2), num_clusters);
		parfor i=1:size(projections, 2)
		    distVec = repmat(projections(:,i), [1 num_clusters]) - C';
		    D(i,:) = sum(distVec.^2);
		end

		percentile95 = prctile(min(D,[],2), 90);
		filteredY = projections(:,min(D,[],2) < percentile95);
		filteredIdx = idx(min(D,[],2) < percentile95);
        filtered_theta = original_theta(min(D,[],2) < percentile95);

		fprintf(1,'Number of Projections Filtered: %d, out of Total: %d', size(projections,2)-size(filteredY,2),size(projections,2));
		clusteredProj = zeros(size(projections,1), num_clusters);

		count = 1;
		% ----------------------------------------------------------------------

		% AVERAGING
		% ----------------------------------------------------------------------
		for ang = 1:num_clusters
		    clusterProjs = filteredY(:,filteredIdx == ang);
		    % If some clusters are obtained
		    if(~(size(clusterProjs,2) == 0))
		        clusteredProj(:,count) = mean(clusterProjs,2);
		        clustered_angles(ang) = mean(filtered_theta(find(filteredIdx == ang)));
		        if count == 1
		            cluster1 = clusterProjs;
		        end
		        if count == 2
		            cluster2 = clusterProjs;
		        end
		        count = count+1;
		    end
		end

		% % Visualization of outliers
		% for k=1:num_clusters
		%     accepted_projections = projections(:, (idx == k) & (min(D,[],2) < percentile95));
		%     rejected_projections = projections(:, (idx == k) & (min(D,[],2) >= percentile95));
		%     rejected_indices = find((idx == k) & (min(D,[],2) >= percentile95));

		%     if size(rejected_projections, 2) > 0 && size(accepted_projections, 2) > 10 && size(rejected_projections, 2) < 3

		%         if sum(ismember(outlier_indices, rejected_indices)) > 0
		%             figure('DefaultAxesFontSize',18)
		%             img_idx = find(ismember(outlier_indices, rejected_indices), 1);
		% %             imshow(allOutlierImages(:, :, img_idx), []);
		%             hold on;
		%             p1 = plot(accepted_projections, 'c');
		%             p2 = plot(rejected_projections, 'r', 'LineWidth', 1);
		%             p3 = plot(clusteredProj(:, k), 'k', 'LineWidth', 2);
		%             hold off;
		%             h = [p1(1);p2(1);p3(1)];
		%             legend(h,'Accepted projections','Rejected Projections','Averaged Projection');
		%         end
		%     end
		% end

		% Trim out the non-existant angles
		clustered_projections = clusteredProj(:,1:count-1);
    end
		
end