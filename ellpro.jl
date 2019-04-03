# Volume properties of an allipsoid
# =================================
# densite defines the mass density;
# ellipsoid is defined by the centre x, radii a, b, c, and orientation R;
# the subroutine outpuss mass, volume, Euler tensor E, and inertia tensor J;
# __________________________________________________________________________
function ellipsoid_properties (a, b, c, R, density)

  volume = (4.0/3.0) * pi * a * b * c;
  mass = density * volume;
  
  J = [mass*(b*b + c*c)/5.0 0.0 0.0;
       0.0 mass*(a*a + c*c)/5.0 0.0;
       0.0 0.0 mass*(b*b + a*a)/5.0];
  
  # note that Inertia = Trace(Euler)*Identity - Euler,
  # hence Euler = 0.5 * Trace (Inertia)*Identity - Inertia,
  # as Trace(Inertia) = 3*Trace(Euler) - Trace(Euler)
  
  E = 0.5*trace(J)*eye(3,3) - J;
  
  E = R*E0*R';

  return mass, volume, E, J
end
