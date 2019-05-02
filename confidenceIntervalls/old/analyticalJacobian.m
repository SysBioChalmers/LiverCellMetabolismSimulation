syms t u f X0 S0 real
anJac = jacobian(S0 + (f/u)*X0*(exp(u*t)-1) , [u, f])
conjSquared = anJac'*anJac

covData = inv(anJac'*anJac)

