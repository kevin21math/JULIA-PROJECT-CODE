using FastGaussQuadrature

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

function galerkin_cross(basis_test, basis_trial, ops::Tuple; a=nothing, b=nothing, nquad=16)
    op_test, op_trial = ops
    if a === nothing || b === nothing
        kv1 = BSplineKit.knots(basis_test); kv2 = BSplineKit.knots(basis_trial)
        a = a === nothing ? min(kv1[1], kv2[1]) : a
        b = b === nothing ? max(kv1[end], kv2[end]) : b
    end
    nb_test = length(basis_test)
    nb_trial = length(basis_trial)
    M = zeros(Float64, nb_test, nb_trial)
    ξ_ref, ω_ref = gausslegendre(nquad)
    mid = 0.5*(a + b)
    half = 0.5*(b - a)
    xq = mid .+ half .* ξ_ref
    wq = half .* ω_ref
    for q in 1:length(xq)
        x = xq[q]; w = wq[q]
        vals_test = [eval_op_scalar(basis_test, i, x, op_test) for i in 1:nb_test]
        vals_trial = [eval_op_scalar(basis_trial, j, x, op_trial) for j in 1:nb_trial]
        M .+= w .* (vals_test * transpose(vals_trial))
    end
    return M
end
