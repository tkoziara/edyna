# Conversion of rotation matrix to quaternion
# ===========================================
# References
# ----------
# [1] http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche52.html
# --------------------------------------------------------------------------
# rotation matrix R is passed as input; quaternion q is returned
#_______________________________________________________________

function rotqua(R)

  q0 = ( R[1,1] + R[2,2] + R[3,3] + 1.0) / 4.0
  q1 = ( R[1,1] - R[2,2] - R[3,3] + 1.0) / 4.0
  q2 = (-R[1,1] + R[2,2] - R[3,3] + 1.0) / 4.0
  q3 = (-R[1,1] - R[2,2] + R[3,3] + 1.0) / 4.0

  if q0 < 0.0 q0 = 0.0 end
  if q1 < 0.0 q1 = 0.0 end
  if q2 < 0.0 q2 = 0.0 end
  if q3 < 0.0 q3 = 0.0 end

  q0 = sqrt(q0)
  q1 = sqrt(q1)
  q2 = sqrt(q2)
  q3 = sqrt(q3)

  if q0 >= q1 && q0 >= q2 && q0 >= q3
      q0 *= 1.0
      q1 *= sign(R[3,2] - R[2,3])
      q2 *= sign(R[1,3] - R[3,1])
      q3 *= sign(R[2,1] - R[1,2])
  elseif q1 >= q0 && q1 >= q2 && q1 >= q3
      q0 *= sign(R[3,2] - R[2,3])
      q1 *= 1.0
      q2 *= sign(R[2,1] + R[1,2])
      q3 *= sign(R[1,3] + R[3,1])
  elseif q2 >= q0 && q2 >= q1 && q2 >= q3
      q0 *= sign(R[1,3] - R[3,1])
      q1 *= sign(R[2,1] + R[1,2])
      q2 *= 1.0
      q3 *= sign(R[3,2] + R[2,3])
  elseif q3 >= q0 && q3 >= q1 && q3 >= q2
      q0 *= sign(R[2,1] - R[1,2])
      q1 *= sign(R[3,1] + R[1,3])
      q2 *= sign(R[3,2] + R[2,3])
      q3 *= 1.0
  else 
      println("rotation_to_quaternion coding error")
  end

  r = sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);

  q0 /= r
  q1 /= r
  q2 /= r
  q3 /= r

  return [q0; q1; q2; q3]

end
