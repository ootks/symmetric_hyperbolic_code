using DynamicPolynomials
using Combinatorics
using CSDP
using SumOfSquares
using PermutationGroups

d = 6
n = d
for i in 1:10
    # Pick d-1 points, r1...rd-1 > 0
    roots = rand(Float64, d-1)
    # roots = [1,1]
    # Normalize them so that their sum is 1.
    roots = roots / sum(roots)

    # Make a polynomial q(t) = (t+1)(t-1/r1)(t-1/r2)...(t-1/rd-1) = c1 + c2t^2+c3t^3+...
    @polyvar t
    uni = prod([t-1/r for r in roots]) * (t+1)
    # Make the polynomial r(t) = q'(t)/t
    reduced = div(differentiate(uni, t), t)

    # Recursively make the polynomial
    # p = 1/2! r(-1) e1^(d-2)e2 + 1/3! r'(-1) e1^(d-3)e3 + ... 
    @polyvar x[1:n]
    function elem(x, k)
        return sum(prod(s) for s in combinations(x,k))/binomial(length(x),k)
    end
    es = [elem(x, k) for k=1:d]
    poly = 0
    for i in 2:d
        poly -= reduced(t=>-1) * es[1]^(d-i) * es[i] / i
        reduced = differentiate(reduced, t) / i
    end

    # Check that the mixed derivative is SOS.
    function wronskian(p, v, w)
        grad = [differentiate(p, x[i]) for i in 1:n]
        diff1 = sum([v[i] * grad[i] for i in 1:n])
        diff2 = sum([w[i] * grad[i] for i in 1:n])
        double_diff = sum([w[i] * differentiate(diff1, x[i]) for i in 1:n])
        return diff1 * diff2 - p * double_diff
    end

    wronsk = wronskian(poly, [1 for i=1:n], [1 for i=1:n])

    solver = optimizer_with_attributes(CSDP.Optimizer)
    model = SOSModel(solver)
    @constraint(model, cref, wronsk >= 0)
    @objective(model, Max, 0)
    JuMP.optimize!(model)
end
