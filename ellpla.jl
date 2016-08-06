# Contact between an ellipsoid and a plane
# ========================================
# References
# ----------
# [1] P. P. Klein, On the Ellipsoid and Plane Intersection
#     Equation, Applied Mathematics, 2012, 3, 1634-1640
# -----------------------------------------------------------------------
# ellipsoid is defined by the centre x, radii a, b, c, and orientation R;
# plane is defined by point and normal;
# _____________________________________________________________
function ellipsoid_plane_contact (x, a, b, c, R, point, normal)
 
  # these are 3-vectors 
  x = vec(x);
  point = vec(point);
  normal = vec(normal);

  # we move the plane into the natural coordinate system of our input ellipsoid,
  # which is centered at x and oriented according to the R rotation
  
  point = R'*(point - x); # so that plane's point, relative to the origin of
                          # the ellipsoid is suitably shifted and rotated
                          # and this needs to be followed by

  normal = R' * normal; # to redefine the normal in the natural coordinate
                        # system of the ellipsoid

  C1 = [1.0/a 0.0 0.0; 0.0 1.0/b 0.0; 0.0 0.0 1.0/c]; # this is matrix D1 in [1]
   
  v = vec(rand(3,1));   # let's generate a random vector first
  r = cross(normal, v); # now r'*normal = 0 as required
  r = r / norm(r);      # lets normalize it (make it's length = 1)
  s = cross(normal, r); # now s'*normal = 0 and s'*r = 0 as required
  q = point;
   
  # now lets check for condition (7) to be fullfilled 
  
  # by doing this transformation we are already satisfying condition (4-6)
  # as rr and ss are orthogonal to each other and to the normal vector
  
  # now lets choose omega and fulfill condition (7)
  
  # here I will define all the dot products needed
  D1 = dot((C1*r),(C1*s));
  D2 = dot((C1*r),(C1*r));
  D3 = dot((C1*s),(C1*s));
  
  if abs(D2-D3) < 1E-10
    omega = 0.25*pi;
  else
    omega = 0.5 * atan((2*D1)/(D2-D3));
  end
  
  rr = cos(omega)* r + sin(omega) * s;  # here we are transforming r and
  ss = -sin(omega) * r + cos(omega) * s; # s to fullfill condition (7) in [1]
  
  r = rr;
  s = ss;
  
  D2 = dot((C1*r),(C1*r));
  D3 = dot((C1*s),(C1*s));
  D4 = dot((C1*q),(C1*r));
  D5 = dot((C1*q),(C1*s));
  D6 = dot((C1*q),(C1*q));
  
  # now lets find the ellipse of the intersection
  d = D6 -((D4)^2/D2)-((D5)^2/D3); # this is just a value
  
  # semi-axes squared radii of the contact ellipse
  Aa = (1-d)/D2;
  Bb = (1-d)/D3;
  
  if Aa > 0. && Bb > 0. # there is contact
    # radii of the contact ellipse
    A  = sqrt(Aa);
    B  = sqrt(Bb);
  
    # now the centre of the contact ellipse
    t0 = -D4/D2;
    u0 = -D5/D3;
    x0 = q + t0*r + u0*s;
  
    # now lets find the contact point (in ellipsoid coordinates)
    alpha = sqrt(1/((x0(1)^2/a^2) + (x0(2)^2/b^2) + (x0(3)^2/c^2)));
    conpnt = alpha * x0;
    
    # and the penetration depth
    depth = dot(normal,  x0 - conpnt);
    
    # now contact point in input coordinates
    conpnt = R1*conpnt + x;

  else # there is no contact
    conpnt = [0.0; 0.0; 0.0];
    depth = -1.0;
    A = 0.0;
    B = 0.0;
  end 

  return conpnt, depth, A, B
end
