function p_d = process_deltas(delta)
    p_d = zeros(size(delta));
    for i=1:size(delta, 2)
        if delta(i) < 0
            p_d(i) = -delta(i) - 110;
        elseif delta(i) > 0
            p_d(i) = delta(i) - 70;
        end
    end
end