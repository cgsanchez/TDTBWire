function generate_wire(N)
    R = Matrix{Float64}(undef, N, 3)
    for i = 0:N-1
        R[i+1, 1] = i * dnn
        R[i+1, 2] = 0.0
        R[i+1, 3] = 0.0
    end
    # Calculate geometrical center of the structure
    Rc = mean(R, dims = 1)
    # Shift the wire to the origin
    R = R .- Rc
    return R
end