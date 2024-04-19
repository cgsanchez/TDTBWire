function generate_H0_from_coordinates(R)
    n_points = size(R, 1)
    H = zeros(n_points, n_points)
    for i = 1:n_points
        for j = i+1:n_points
            distance = norm(R[i, :] .- R[j, :])
            if distance <= dnn
                H[i, j] = H[j, i] = β
            end
        end
    end
    return H
end