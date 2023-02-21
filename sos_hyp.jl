using SumOfSquares
using DynamicPolynomials
using CSDP
import Polynomials
using Combinatorics
using LinearAlgebra

using Random

n = 5
d = 5

println("Starting!")
@polyvar x[1:n]
function elem(x, k)
    return sum(prod(s) for s in combinations(x,k))
end
es = [elem(x, k) for k=1:d]
# gammas = [1,0,23/10,33/10,68/15,6]
es_coeffs = [0.,0.,7.,-220.,4500.]
p = sum(es_coeffs[i]*es[1]^(d-i)*es[i] for i=1:d)
v = [6, 1.,1.,1.,1.]

function normalize(p, v)
    @polyvar t
    univar = p([x[i] => t+v[i] for i in 1:n]...)
    coeffs = [DynamicPolynomials.coefficient(t) for t in terms(univar)]
    univar2 = Polynomials.Polynomial(reverse(coeffs))
    roots = (Polynomials.roots(univar2))
    if !all([isreal(root) for root in roots])
        error("Not real rooted")
    end
    root = max(roots...)
    return [root + v[i] for i in 1:n]
end

function wronskian(p, v, w)
    grad = [differentiate(p, x[i]) for i in 1:n]
    diff1 = sum([v[i] * grad[i] for i in 1:n])
    diff2 = sum([w[i] * grad[i] for i in 1:n])
    double_diff = sum([w[i] * differentiate(diff1, x[i]) for i in 1:n])
    return diff1 * diff2 - p * double_diff
end
wronsk = wronskian(p, [1,1,1,1,1.], [1,1,1,1,1.])
solver = optimizer_with_attributes(CSDP.Optimizer)
model = SOSModel(solver)
@constraint(model, cref, wronsk >= 0)
@objective(model, Max, 0)
JuMP.optimize!(model)

print(sos_decomposition(cref, 1e-4))

#wronsk = wronskian(p, v, [1.,1,1,1,1])
#solver = optimizer_with_attributes(CSDP.Optimizer)
#model = SOSModel(solver)
#@constraint(model, cref, wronsk >= 0)
#@objective(model, Max, 0)
#JuMP.optimize!(model)

function rational(x, k) 
    if k == 1
        return floor(Int, x)
    end
    int_part = floor(Int, x)
    return int_part + 1//rational(1/(x-int_part), k-1)
end
rat = x -> convert(Rational{BigInt}, rational(x, 10))
x = rat.(convert(Array,moment_matrix(cref).Q))
basis = moment_matrix(cref).basis.monomials

function matrixize(f, basis)
    taken = Set()
    matrix = zeros(length(basis),length(basis))
    term_list = terms(f)
    term_dict = Dict(monomial(term) => DynamicPolynomials.coefficient(term) for term in term_list)
    for (i, mono1) in enumerate(basis)
        for (j, mono2) in enumerate(basis)
            prod = mono1 * mono2
            if in(prod, taken)
                continue
            end
            push!(taken, prod)
            if in(prod, keys(term_dict))
                matrix[i,j] = term_dict[prod]
            end
        end
    end
    return matrix
end
y = matrixize(wronsk, basis)
println(y)
println(tr(x * y))
#
#function check_sos(p, v, w)
#    wronsk = wronskian(p, v, w)
#    println("Writing Wronskian")
#    open("wronskian.txt", "w") do file
#        write(file, string(wronsk))
#    end
#    solver = optimizer_with_attributes(CSDP.Optimizer)
#    model = SOSModel(solver)
#    @constraint(model, cref, wronsk >= 0)
#    @objective(model, Max, 0)
#    JuMP.optimize!(model)
#    println("Writing Moments")
#    open("moments.txt", "w") do file
#        write(file, string(moments(cref)))
#    end
#    println("Writing Moments matrix")
#    open("moments_matrix.txt", "w") do file
#        write(file, string(moment_matrix(cref)))
#    end
#    println(eigen(convert(Array,moment_matrix(cref).Q)).values)
#    println("Writing dual")
#    open("dual.txt", "w") do file
#        write(file, string(dual(cref)))
#    end
#    return Int(primal_status(model)) != 1
#end
##check_sos(es[5], [1,1,1,1,1.],[1,1,1,1,1.])
#
#v = [6, 1.,1.,1.,1.]
#check_sos(p, [1.,1,1,1,1], v)
##iters = 10
##for a in (-iters):iters
##    for b in (-iters):a
##        for c in (-iters):b
##            for d in (-iters):c
##                v = [1, a/iters, b/iters, c/iters, d/iters]
##                try 
##                    v = normalize(p,v)
##                catch e
##                    println("not real rooted")
##                    println(v)
##                    continue
##                end
##                w = [1.,1,1,1,1]
##                if check_sos(p, v, w)
##                    println("Not SOS!")
##                    println("v: ", v)
##                    println("Wronskian: ", wronsk)
##                    break
##                end
##            end
##		end
##    end
##end
