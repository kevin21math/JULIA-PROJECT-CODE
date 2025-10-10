using FastGaussQuadrature
using BSplineKit


function galerkin_cross(basis_test::BSplineBasis, basis_trial::BSplineBasis,
                        ops::Tuple{BSplineKit.Derivative, BSplineKit.Derivative};
                        nquad::Int = 8)
    Dtest, Dtrial = ops
    knots1 = unique(BSplineKit.knots(basis_test))
    knots2 = unique(BSplineKit.knots(basis_trial))
    a = min(knots1[1], knots2[1])
    b = max(knots1[end], knots2[end])

    knots = sort(unique(vcat(knots1, knots2)))
    nb_test = length(basis_test)
    nb_trial = length(basis_trial)
    M = zeros(Float64, nb_test, nb_trial)
    ξ_ref, ω_ref = gausslegendre(nquad)

    for e = 1:(length(knots)-1)
        xL, xR = knots[e], knots[e+1]
        if xR <= xL
            continue
        end
        mid = 0.5*(xL + xR)
        half = 0.5*(xR - xL)
        xq = mid .+ half .* ξ_ref
        wq = half .* ω_ref
        nq = length(xq)
        for q in 1:nq
            x = xq[q]; w = wq[q]
            vals_test  = [evaluate(basis_test, i, x, Dtest)  for i in 1:nb_test]
            vals_trial = [evaluate(basis_trial, j, x, Dtrial) for j in 1:nb_trial]
            M .+= w .* (vals_test * transpose(vals_trial))
        end
    end
    return M
end
