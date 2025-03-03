function lhs = latin_hypercube_sampling(dim, size, num_samples)
    lhs = zeros(dim, num_samples);
    for i = 1:num_samples
        lhs(:, i) = (lhsdesign(dim, 1) - 0.5) * size;
    end
end