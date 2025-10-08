using FastGaussQuadrature, SparseArrays

function eval_op_scalar(basis, j::Integer, x::Real, op)
    if op === :D0 || op === nothing
        return Float64(BSplineKit.BSplines.evaluate(basis, j, x))
    end
    try
        return Float64(BSplineKit.BSplines.evaluate(basis, j, x, op))
    catch _
        h = 1e-8 * max(1.0, abs(x))
        vp = Float64(BSplineKit.BSplines.evaluate(basis, j, x + h))
        vm = Float64(BSplineKit.BSplines.evaluate(basis, j, x - h))
        return (vp - vm) / (2h)
    end
end

function galerkin_cross_elementwise(basis_test, basis_trial, ops::Tuple; nquad=8)
    op_test, op_trial = ops
    kv1 = BSplineKit.knots(basis_test)
    kv2 = BSplineKit.knots(basis_trial)
    breaks = unique(sort(vcat(kv1, kv2)))
    elems = [(breaks[i], breaks[i+1]) for i in 1:(length(breaks)-1) if breaks[i+1] > breaks[i] + 1e-14]
    nb_test = length(basis_test)
    nb_trial = length(basis_trial)
    I = Int[]
    J = Int[]
    V = Float64[]
    ξ_ref, ω_ref = gausslegendre(nquad)
    for (a,b) in elems
        mid = 0.5*(a + b)
        half = 0.5*(b - a)
        xq = mid .+ half .* ξ_ref
        wq = half .* ω_ref
        idx_test = [i for i in 1:nb_test if any(abs(eval_op_scalar(basis_test,i,x,op_test)) > 0 for x in xq)]
        idx_trial = [j for j in 1:nb_trial if any(abs(eval_op_scalar(basis_trial,j,x,op_trial)) > 0 for x in xq)]
        if isempty(idx_test) || isempty(idx_trial)
            continue
        end
        Mloc = zeros(Float64, length(idx_test), length(idx_trial))
        for q in 1:length(xq)
            x = xq[q]; w = wq[q]
            vals_t = [eval_op_scalar(basis_test, i, x, op_test) for i in idx_test]
            vals_r = [eval_op_scalar(basis_trial, j, x, op_trial) for j in idx_trial]
            Mloc .+= w .* (vals_t * transpose(vals_r))
        end
        for ii in 1:length(idx_test), jj in 1:length(idx_trial)
            val = Mloc[ii,jj]
            if abs(val) > 1e-14
                push!(I, idx_test[ii])
                push!(J, idx_trial[jj])
                push!(V, val)
            end
        end
    end
    return sparse(I, J, V, nb_test, nb_trial)
end
