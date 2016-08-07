# Exponential map via Rodriguez formula
# =====================================
function expmap(O)
  O = vec(O);
  dp = dot(O,O);
  ln = sqrt(dp);
  W = skew(O);
  if ln < 1E-15
    R = eye(3,3);
  else
    R = eye(3,3) + (sin(ln)/ln)*W + ((1-cos(ln))/dp)*W*W;
  end
  return R;
end

# Skew-symmetric matrix out of 3-vector
# =====================================
function skew(O)
  [0.0 -O[3] O[2];
   O[3] 0.0 -O[1];
  -O[2] O[1] 0.0];
end 
