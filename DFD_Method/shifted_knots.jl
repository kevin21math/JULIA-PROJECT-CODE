
function shifted_knots(a::Real, b::Real, N::Integer, p::Integer, αk::Real)
    @assert a < b 
    @assert N >= p
    @assert p >= 1

    α = clamp(float(αk), -1.0, 1.0)
    order = p
    nk = N + order                
    n_internal = nk - 2*order    
    internal = zeros(Float64, n_internal)
   
    if n_internal > 0
        for i = 1:n_internal
           internal[i] = a + ((b-a)/2)*(1-(N+p-2*(p+i-1))/(N-1))
        end
        #internal = collect(range(a, b, length = n_internal + 2))[2:end-1]  
    else
        internal = Float64[]
    end

    
    kv = zeros(Float64, nk)
    for i in 1:order
        kv[i] = a
    end
    for i in 1:n_internal
        kv[order + i] = internal[i]
    end
    for i in (order + length(internal) + 1):nk
        kv[i] = b
    end

    if n_internal > 0
        for idx in (order+1):(order + n_internal)
            local_i = idx - order             
            curr = internal[local_i]
            right = (local_i < n_internal) ? internal[local_i + 1] : b
            left  = (local_i > 1)            ? internal[local_i - 1] : a
            if α >= 0
                kv[idx] = (1 - α) * curr + α * right
            else
                kv[idx] = (1 - α) * curr + α * left
            end
        end
    end

    return kv
end
