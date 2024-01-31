# Application of the function H to a Mueller matrix M
# The function is taken from the paper 
# Jose J. Gil, "On optimal filtering of measured Mueller matrices", arXiv:1605.04704


function H(M)
  return [M[1,1] + M[1,2] + M[2,1] + M[2,2]         M[1,3] + M[2,3] + im*(M[1,4] + M[2,4])   M[3,1] + M[3,2] - im*(M[4,1] + M[4,2])   M[3,3] + M[4,4] + im*(M[3,4] - M[4,3])
          M[1,3] + M[2,3] - im*(M[1,4] + M[2,4])    M[1,1] - M[1,2] + M[2,1] - M[2,2]        M[3,3] - M[4,4] - im*(M[3,4] + M[4,3])   M[3,1] - M[3,2] - im*(M[4,1] - M[4,2])
          M[3,1] + M[3,2] + im*(M[4,1] + M[4,2])    M[3,3] - M[4,4] + im*(M[3,4] + M[4,3])   M[1,1] + M[1,2] - M[2,1] - M[2,2]        M[1,3] - M[2,3] + im*(M[1,4] - M[2,4])
          M[3,3] + M[4,4] - im*(M[3,4] - M[4,3])    M[3,1] - M[3,2] + im*(M[4,1] - M[4,2])   M[1,3] - M[2,3] - im*(M[1,4] - M[2,4])   M[1,1] - M[1,2] - M[2,1] + M[2,2]
          ]

end
