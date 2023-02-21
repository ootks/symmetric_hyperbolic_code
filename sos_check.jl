using DynamicPolynomials
using TypedPolynomials
using SumOfSquares
import CSDP


DynamicPolynomials.@polyvar x1 x2 x3 x4

r = 4*x1^2*x2^2+2*x1^2*x2*x3+2*x1*x2^2*x3+8*x1^2*x3^2+2*x1*x2*x3^2+6*x2^2*x3^2+2*x1^2*x2*x4+2*x1*x2^2*x4+6*x1^2*x3*x4+4*x2^2*x3*x4+6*x1*x3^2*x4+4*x2*x3^2*x4+8*x1^2*x4^2+2*x1*x2*x4^2+6*x2^2*x4^2+6*x1*x3*x4^2+4*x2*x3*x4^2+10*x3^2*x4^2




model = SOSModel(CSDP.Optimizer)
@constraint(model, cref, r >= 0)
optimize!(model)
sos_dec = sos_decomposition(cref, 1e-4)