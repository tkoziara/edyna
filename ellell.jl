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

  A = R1*[1.0/a1^2 0.0 0.0; 0.0 1.0/b1^2 0.0; 0.0 0.0 1.0/c1^2]*R1';
  B = R2*[1.0/a2^2 0.0 0.0; 0.0 1.0/b2^2 0.0; 0.0 0.0 1.0/c2^2]*R2';

  return conpnt, depth, A, B
end
