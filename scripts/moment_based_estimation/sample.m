function keep = sample(P, numkeep)
    keep = randperm(size(P,2), numkeep);



    % 
    %     % OR
    %     % Sample from a skewed distribution
    %     weights = [0.2 0.3 0.12 0.03 0.35]; 
    %     partitionsize = length(thetas)/length(weights);
    %     
    %     keep = [];
    %     for i = 1:length(weights)
    %         start = partitionsize * (i - 1);
    %         numsample = round(numkeep * weights(i));
    %         keep_partition = start + randperm(partitionsize, numsample);
    %         keep = [keep  keep_partition];
    %     end
    
end