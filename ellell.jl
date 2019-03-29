include("ellpla.jl")
# Contact between two ellipsoids
# ==============================
# References
# ----------
# [1] J. W. Perram, ￼Statistical Mechanics of Hard Ellipsoids. I.
#     Overlap Algorithm and the Contact Function,
#     Journal of Computational Physics, 58, 409-416 (1985)
# -----------------------------------------------------------------------
# ellipsoids are defined by centre x, radii a, b, c, and orientation R;
# -----------------------------------------------------------------------
# returned: (conpnt, depth, A, B), where:
# conpnt - contact point (or [0.0, 0.0, 0.0] if no contact)
# depth - positive penetration depth (or -1.0 if no contact)
# A, B - approximate contact ellipse radii (or 0.0, 0.0 if no contact)
# __________________________________________________________________________
function ellipsoid_ellipsoid_contact(x1, a1, b1, c1, R1, x2, a2, b2, c2, R2)

  x1 = vec(x1);
  x2 = vec(x2);

  IA = R1*[1.0/(a1*a1) 0.0 0.0; 0.0 1.0/(b1*b1) 0.0; 0.0 0.0 1.0/(c1*c1)]*R1';
  IB = R2*[1.0/(a2*a2) 0.0 0.0; 0.0 1.0/(b2*b2) 0.0; 0.0 0.0 1.0/(c2*c2)]*R2';

  A = R1*[a1*a1 0.0 0.0; 0.0 b1*b1 0.0; 0.0 0.0 c1*c1]*R1';
  B = R2*[a2*a2 0.0 0.0; 0.0 b2*b2 0.0; 0.0 0.0 c2*c2]*R2';

  x21 = x2 - x1;

  C = Matrix{Float64}(undef, 3, 3)
  x = Vector{Float64}(undef, 3)
  lambda = 0.5
  error = 1.0
  iter = 0

  while error > 1E-10 && iter < 100

    C = inv(lambda*B + (1.0-lambda)*A)
    x = x2 - lambda*B*C*x21
    #x = x1 + (1.0-lambda)*A*C*x21 # alternative formula
    #x0 = (lambda*IA + (1.0-lambda)*IB)\(lambda*IA*x1 + (1.0-lambda)*IB*x2) # hand calculation
    #println("x0 = ", x0, ", x = ", x, "x - x0 = ", x-x0) # comparison with paper formula
    FA = (x-x1)'*IA*(x-x1)
    FB = (x-x2)'*IB*(x-x2)
    dFdl = FA - FB
    #ddFdll = -2.0*x21'*C*(lambda*IA + (1.0-lambda)*IB)*C*x21 # paper formula, poorer rate
    #ddFdll = -2.0*x21'*C*(lambda*A + (1.0-lambda)*B)*C*x21 # why this is better?
    #ddFdll = -2.0*x21'*C*(lambda*B + (1.0-lambda)*A)*C*x21 # why this is so much better?
    ddFdll = -2.0*x21'*C*x21 # (same as above) why this is so much better?
    dlambda = -dFdl/ddFdll
    lambda += dlambda

    error = abs(dlambda)
    iter += 1

    println("at iter ", iter, " error ", error, ", and x = ", x)
  end

  F = lambda*(1.0-lambda)*x21'*C*x21

  if F < 1.0
    v = lambda*IA*(x-x1) - (1.0-lambda)*IB*(x-x2)
    conpnt = x
    normal = v/norm(v)
    (y1, d1, A1, B1) = ellipsoid_plane_contact(x1, a1, b1, c1, R1, conpnt, -normal)
    (y2, d2, A2, B2) = ellipsoid_plane_contact(x2, a2, b2, c2, R2, conpnt, normal)
    println(y1, " ", d1, " ", A1, " ", B1)
    println(y2, " ", d2, " ", A2, " ", B2)
    depth = d1+d2
    A = min(A1, A2)
    B = min(B1, B2)
  else
    conpnt = [0.0; 0.0; 0.0];
    normal = [0.0; 0.0; 0.0];
    depth = -1.0;
    A = 0.0;
    B = 0.0;
  end

  return conpnt, normal, depth, A, B
end
