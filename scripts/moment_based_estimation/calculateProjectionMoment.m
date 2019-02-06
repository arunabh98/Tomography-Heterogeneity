function M = calculateProjectionMoment(P, svector, k)
    % For a vector P
    % Calculate the k^th moment of the data P.
    % The indices in P are to be taken from svector.
    
    M = 0;
    for i = 1:length(P)
        addition = P(i) * (svector(i)^k);
        M = M + addition;
    end
end
