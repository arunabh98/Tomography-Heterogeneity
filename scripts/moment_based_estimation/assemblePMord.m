function PMord = ...
    assemblePMord(Pgiven, kmax, svector, Ord, numkeep)
    M = ones(size(Pgiven,2), kmax+1); % each column represents kth moments for one k, for all directions
    for k = 0:kmax
        for i = 1:size(Pgiven,2)
            M(i, k+1) = calculateProjectionMoment(Pgiven(:, i), svector, k);
        end
    end
    PMord = reshape(M(:,2:(Ord+1)), Ord * numkeep, 1);
end