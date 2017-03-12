include("ellpla.jl")
# -------------------------------------------
# Contact between an ellipsoid and a triangle
# ===========================================
# References
# -----------
# [1] http://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates
# -----------------------------------------------------------------------------------------------------------------
# ellipsoid is defined by centre x, radii a, b, c, and orientation R;
# triangle is identified by points p, q, r;
# ----------------------------------------------
# returned (conpnt, depth, A, B, normal), where:
# conpnt - contact point (or [0.0, 0.0, 0.0] if no contact)
# normal - contact normal (or [0.0, 0.0, 0.0] if no contact)
# depth - positive penetration depth (or -1.0 if no contact)
# A, B - contact ellipse radii (or 0.0, 0.0 if no contact)
# ____________________________________________________________
function ellipsoid_triangle_contact(x, a, b, c, R, p, q, r)
  e0 = q-p
  e1 = r-p
  normal = cross(-e1, e0)
  normal = normal/norm(normal)
  (conpnt, depth, A, B) = ellipsoid_plane_contact(x, a, b, c, R, p, normal);
  if depth == -1.0
    normal = [0.0 0.0 0.0];
  else
    vv = conpnt-p;
    aa = dot(e0,e0)
    bb = dot(e0,e1)
    cc = dot(e1,e1)
    dd = dot(e0,vv)
    ee = dot(e1,vv)
    mm = aa*cc-bb*bb
    ss = (cc*dd-bb*ee)/mm
    tt = (aa*ee-bb*dd)/mm
    if !(ss >= 0.0 && tt >= 0.0 && ss + tt <= 1.0)
      conpnt = [0.0; 0.0; 0.0];
      normal = [0.0; 0.0; 0.0];
      depth = -1.0;
      A = 0.0;
      B = 0.0;
    end
  end
  return (conpnt, normal, depth, A, B)
end
