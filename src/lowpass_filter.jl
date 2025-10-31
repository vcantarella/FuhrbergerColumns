function lowpass_filter(x::Vector{Float64}, α::Float64, )
    y = similar(x)
    y[1] = x[1]
    for n in 2:length(x)
        y[n] = α  * y[n - 1] + (1 - α) * x[n]
    end
    return y
end

function predict!(x, A, P, Q)
    # Predict the next state
    x .= A * x
    # Predict the next covariance
    P .= A * P * A' + Q
end

function update!(x, P, z, H, R)
    # Compute the Kalman gain
    K = P * H' \ (H * P * H' + R)
    # Update the state estimate
    x .= x + K * (z - H * x)
    # Update the covariance estimate
    P .= (I - K * H) * P
end