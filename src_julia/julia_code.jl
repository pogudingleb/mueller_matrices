using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Folds

# Size
sz = 600*700

# Download the input data into a matrix of arrays
println()
println("Benchmarking the performance of verifying positive definite and positive semidefinite properties using Choletsky, leading principal minors and eigenvalues")
println()
println("Time to download the input data into a matrix of arrays:")
@btime global MM = [reshape(readdlm("data/MM_"*string(i)*string(j)*".dat"), 1, sz) for i in 1:4, j in 1:4];
println()


# Convert the matrix of arrays into an array of matrices

println("Time to convert the matrix of arrays into an array of matrices:")
@btime global MMlist = [[MM[i,j][s] for i in 1:4, j in 1:4] for s in 1:sz];
println()

# The function H from Gil's paper
function H(M)
  return [M[1,1] + M[1,2] + M[2,1] + M[2,2]         M[1,3] + M[2,3] + im*(M[1,4] + M[2,4])   M[3,1] + M[3,2] - im*(M[4,1] + M[4,2])   M[3,3] + M[4,4] + im*(M[3,4] - M[4,3])
          M[1,3] + M[2,3] - im*(M[1,4] + M[2,4])    M[1,1] - M[1,2] + M[2,1] - M[2,2]        M[3,3] - M[4,4] - im*(M[3,4] + M[4,3])   M[3,1] - M[3,2] - im*(M[4,1] - M[4,2])
          M[3,1] + M[3,2] + im*(M[4,1] + M[4,2])    M[3,3] - M[4,4] + im*(M[3,4] + M[4,3])   M[1,1] + M[1,2] - M[2,1] - M[2,2]        M[1,3] - M[2,3] + im*(M[1,4] - M[2,4])
          M[3,3] + M[4,4] - im*(M[3,4] - M[4,3])    M[3,1] - M[3,2] + im*(M[4,1] - M[4,2])   M[1,3] - M[2,3] - im*(M[1,4] - M[2,4])   M[1,1] - M[1,2] - M[2,1] + M[2,2]
          ]

end

# Functions that are coefficients of the characteristic polynomial

function C1(M)
  return M[1,1]+M[2,2]+M[3,3]+M[4,4]
end

function C2(M)
  return M[1,1]*M[2,2]+M[1,1]*M[3,3]+M[1,1]*M[4,4]-M[1,2]*M[2,1]-M[1,3]*M[3,1]-M[1,4]*M[4,1]+M[2,2]*M[3,3]+M[2,2]*M[4,4]-M[2,3]*M[3,2]-M[2,4]*M[4,2]+M[3,3]*M[4,4]-M[3,4]*M[4,3]
end

function C3(M)
  return M[1,1]*M[2,2]*M[3,3]+M[1,1]*M[2,2]*M[4,4]-M[1,1]*M[2,3]*M[3,2]-M[1,1]*M[2,4]*M[4,2]+M[1,1]*M[3,3]*M[4,4]-M[1,1]*M[3,4]*M[4,3]-M[1,2]*M[2,1]*M[3,3]-M[1,2]*M[2,1]*M[4,4]+M[1,2]*M[2,3]*M[3,1]+M[1,2]*M[2,4]*M[4,1]+M[1,3]*M[2,1]*M[3,2]-M[1,3]*M[2,2]*M[3,1]-M[1,3]*M[3,1]*M[4,4]+M[1,3]*M[3,4]*M[4,1]+M[1,4]*M[2,1]*M[4,2]-M[1,4]*M[2,2]*M[4,1]+M[1,4]*M[3,1]*M[4,3]-M[1,4]*M[3,3]*M[4,1]+M[2,2]*M[3,3]*M[4,4]-M[2,2]*M[3,4]*M[4,3]-M[2,3]*M[3,2]*M[4,4]+M[2,3]*M[3,4]*M[4,2]+M[2,4]*M[3,2]*M[4,3]-M[2,4]*M[3,3]*M[4,2]
end

function C4(M)
  return M[1,1]*M[2,2]*M[3,3]*M[4,4]-M[1,1]*M[2,2]*M[3,4]*M[4,3]-M[1,1]*M[2,3]*M[3,2]*M[4,4]+M[1,1]*M[2,3]*M[3,4]*M[4,2]+M[1,1]*M[2,4]*M[3,2]*M[4,3]-M[1,1]*M[2,4]*M[3,3]*M[4,2]-M[1,2]*M[2,1]*M[3,3]*M[4,4]+M[1,2]*M[2,1]*M[3,4]*M[4,3]+M[1,2]*M[2,3]*M[3,1]*M[4,4]-M[1,2]*M[2,3]*M[3,4]*M[4,1]-M[1,2]*M[2,4]*M[3,1]*M[4,3]+M[1,2]*M[2,4]*M[3,3]*M[4,1]+M[1,3]*M[2,1]*M[3,2]*M[4,4]-M[1,3]*M[2,1]*M[3,4]*M[4,2]-M[1,3]*M[2,2]*M[3,1]*M[4,4]+M[1,3]*M[2,2]*M[3,4]*M[4,1]+M[1,3]*M[2,4]*M[3,1]*M[4,2]-M[1,3]*M[2,4]*M[3,2]*M[4,1]-M[1,4]*M[2,1]*M[3,2]*M[4,3]+M[1,4]*M[2,1]*M[3,3]*M[4,2]+M[1,4]*M[2,2]*M[3,1]*M[4,3]-M[1,4]*M[2,2]*M[3,3]*M[4,1]-M[1,4]*M[2,3]*M[3,1]*M[4,2]+M[1,4]*M[2,3]*M[3,2]*M[4,1]
end

# Function to check positive semi-definite using the coefficients of the characteristic polynomial

function PSD(M)
  return ((real(C1(M)) >= 0) & (real(C2(M)) >= 0) & (real(C3(M)) >= 0) & (real(C4(M)) >= 0))
end


function C1_noH(M)
  return 4*M[1,1]
end

function C2_noH(M)
  return 6*M[1,1]^2-2*M[1,2]^2-2*M[1,3]^2-2*M[1,4]^2-2*M[2,1]^2-2*M[2,2]^2-2*M[2,3]^2-2*M[2,4]^2-2*M[3,1]^2-2*M[3,2]^2-2*M[3,3]^2-2*M[3,4]^2-2*M[4,1]^2-2*M[4,2]^2-2*M[4,3]^2-2*M[4,4]^2
end

function C3_noH(M)
  return 4*M[1,1]^3+(-4*M[1,2]^2-4*M[1,3]^2-4*M[1,4]^2-4*M[2,1]^2-4*M[2,2]^2-4*M[2,3]^2-4*M[2,4]^2-4*M[3,1]^2-4*M[3,2]^2-4*M[3,3]^2-4*M[3,4]^2-4*M[4,1]^2-4*M[4,2]^2-4*M[4,3]^2-4*M[4,4]^2)*M[1,1]+(8*M[2,1]*M[2,2]+8*M[3,1]*M[3,2]+8*M[4,1]*M[4,2])*M[1,2]+(8*M[2,1]*M[2,3]+8*M[3,1]*M[3,3]+8*M[4,1]*M[4,3])*M[1,3]+(8*M[2,1]*M[2,4]+8*M[3,1]*M[3,4]+8*M[4,1]*M[4,4])*M[1,4]+(8*M[3,3]*M[4,4]-8*M[3,4]*M[4,3])*M[2,2]+(-8*M[3,2]*M[4,4]+8*M[3,4]*M[4,2])*M[2,3]+8*M[2,4]*(M[3,2]*M[4,3]-M[3,3]*M[4,2])
end

function C4_noH(M)
  return M[1,1]^4+M[1,2]^4+M[1,3]^4+M[2,1]^4+M[2,2]^4+M[2,3]^4+M[3,1]^4+M[3,2]^4+M[3,3]^4+M[1,4]^4+M[2,4]^4+M[3,4]^4+8*M[3,3]*M[3,4]*M[4,3]*M[4,4]-8*M[4,1]*(M[3,2]*M[4,2]+M[3,3]*M[4,3]+M[3,4]*M[4,4])*M[3,1]+8*M[4,2]*(M[3,3]*M[4,3]+M[3,4]*M[4,4])*M[3,2]+8*M[2,4]*(M[3,3]*M[3,4]+M[4,3]*M[4,4])*M[2,3]+(M[4,1]^2-M[4,2]^2-M[4,3]^2-M[4,4]^2)^2+(2*M[4,1]^2-2*M[4,2]^2-2*M[4,3]^2+2*M[4,4]^2)*M[3,4]^2+(2*M[3,4]^2+2*M[4,1]^2-2*M[4,2]^2+2*M[4,3]^2-2*M[4,4]^2)*M[3,3]^2+(2*M[3,3]^2+2*M[3,4]^2+2*M[4,1]^2+2*M[4,2]^2-2*M[4,3]^2-2*M[4,4]^2)*M[3,2]^2+(-2*M[3,2]^2-2*M[3,3]^2-2*M[3,4]^2+2*M[4,1]^2+2*M[4,2]^2+2*M[4,3]^2+2*M[4,4]^2)*M[3,1]^2+(2*M[3,1]^2-2*M[3,2]^2-2*M[3,3]^2+2*M[3,4]^2+2*M[4,1]^2-2*M[4,2]^2-2*M[4,3]^2+2*M[4,4]^2)*M[2,4]^2+(2*M[2,4]^2+2*M[3,1]^2-2*M[3,2]^2+2*M[3,3]^2-2*M[3,4]^2+2*M[4,1]^2-2*M[4,2]^2+2*M[4,3]^2-2*M[4,4]^2)*M[2,3]^2+((8*M[3,2]*M[3,3]+8*M[4,2]*M[4,3])*M[2,3]+8*M[2,4]*(M[3,2]*M[3,4]+M[4,2]*M[4,4]))*M[2,2]+(2*M[2,3]^2+2*M[2,4]^2+2*M[3,1]^2+2*M[3,2]^2-2*M[3,3]^2-2*M[3,4]^2+2*M[4,1]^2+2*M[4,2]^2-2*M[4,3]^2-2*M[4,4]^2)*M[2,2]^2+((-8*M[3,1]*M[3,2]-8*M[4,1]*M[4,2])*M[2,2]+(-8*M[3,1]*M[3,3]-8*M[4,1]*M[4,3])*M[2,3]-8*M[2,4]*(M[3,1]*M[3,4]+M[4,1]*M[4,4]))*M[2,1]+(-2*M[2,2]^2-2*M[2,3]^2-2*M[2,4]^2+2*M[3,1]^2+2*M[3,2]^2+2*M[3,3]^2+2*M[3,4]^2+2*M[4,1]^2+2*M[4,2]^2+2*M[4,3]^2+2*M[4,4]^2)*M[2,1]^2+((-8*M[3,2]*M[4,3]+8*M[3,3]*M[4,2])*M[2,1]+(8*M[3,1]*M[4,3]-8*M[3,3]*M[4,1])*M[2,2]+8*M[2,3]*(-M[3,1]*M[4,2]+M[3,2]*M[4,1]))*M[1,4]+(-2*M[2,1]^2+2*M[2,2]^2+2*M[2,3]^2-2*M[2,4]^2-2*M[3,1]^2+2*M[3,2]^2+2*M[3,3]^2-2*M[3,4]^2-2*M[4,1]^2+2*M[4,2]^2+2*M[4,3]^2-2*M[4,4]^2)*M[1,4]^2+((-8*M[2,3]*M[2,4]-8*M[3,3]*M[3,4]-8*M[4,3]*M[4,4])*M[1,4]+(8*M[3,2]*M[4,4]-8*M[3,4]*M[4,2])*M[2,1]+(-8*M[3,1]*M[4,4]+8*M[3,4]*M[4,1])*M[2,2]-8*M[2,4]*(-M[3,1]*M[4,2]+M[3,2]*M[4,1]))*M[1,3]+(2*M[1,4]^2-2*M[2,1]^2+2*M[2,2]^2-2*M[2,3]^2+2*M[2,4]^2-2*M[3,1]^2+2*M[3,2]^2-2*M[3,3]^2+2*M[3,4]^2-2*M[4,1]^2+2*M[4,2]^2-2*M[4,3]^2+2*M[4,4]^2)*M[1,3]^2+((-8*M[2,2]*M[2,3]-8*M[3,2]*M[3,3]-8*M[4,2]*M[4,3])*M[1,3]+(-8*M[2,2]*M[2,4]-8*M[3,2]*M[3,4]-8*M[4,2]*M[4,4])*M[1,4]+(-8*M[3,3]*M[4,4]+8*M[3,4]*M[4,3])*M[2,1]+(8*M[3,1]*M[4,4]-8*M[3,4]*M[4,1])*M[2,3]+8*M[2,4]*(-M[3,1]*M[4,3]+M[3,3]*M[4,1]))*M[1,2]+(2*M[1,3]^2+2*M[1,4]^2-2*M[2,1]^2-2*M[2,2]^2+2*M[2,3]^2+2*M[2,4]^2-2*M[3,1]^2-2*M[3,2]^2+2*M[3,3]^2+2*M[3,4]^2-2*M[4,1]^2-2*M[4,2]^2+2*M[4,3]^2+2*M[4,4]^2)*M[1,2]^2+((8*M[2,1]*M[2,2]+8*M[3,1]*M[3,2]+8*M[4,1]*M[4,2])*M[1,2]+(8*M[2,1]*M[2,3]+8*M[3,1]*M[3,3]+8*M[4,1]*M[4,3])*M[1,3]+(8*M[2,1]*M[2,4]+8*M[3,1]*M[3,4]+8*M[4,1]*M[4,4])*M[1,4]+(8*M[3,3]*M[4,4]-8*M[3,4]*M[4,3])*M[2,2]+(-8*M[3,2]*M[4,4]+8*M[3,4]*M[4,2])*M[2,3]+8*M[2,4]*(M[3,2]*M[4,3]-M[3,3]*M[4,2]))*M[1,1]+(-2*M[1,2]^2-2*M[1,3]^2-2*M[1,4]^2-2*M[2,1]^2-2*M[2,2]^2-2*M[2,3]^2-2*M[2,4]^2-2*M[3,1]^2-2*M[3,2]^2-2*M[3,3]^2-2*M[3,4]^2-2*M[4,1]^2-2*M[4,2]^2-2*M[4,3]^2-2*M[4,4]^2)*M[1,1]^2
end

function PSD_noH(M)
  return ((C1_noH(M) >= 0) & (C2_noH(M) >= 0) & (C3_noH(M) >= 0) & (C4_noH(M) >= 0))
end


# Determinants used in Sylvester's criterion

function D1(M)
  return M[1,1]
end

function D2(M)
  return M[1,1]*M[2,2]-M[1,2]*M[2,1]
end

function D3(M)
  return M[1,1]*M[2,2]*M[3,3]-M[1,1]*M[2,3]*M[3,2]-M[1,2]*M[2,1]*M[3,3]+M[1,2]*M[2,3]*M[3,1]+M[1,3]*M[2,1]*M[3,2]-M[1,3]*M[2,2]*M[3,1]
end

function D4(M)
  return M[1,1]*M[2,2]*M[3,3]*M[4,4]-M[1,1]*M[2,2]*M[3,4]*M[4,3]-M[1,1]*M[2,3]*M[3,2]*M[4,4]+M[1,1]*M[2,3]*M[3,4]*M[4,2]+M[1,1]*M[2,4]*M[3,2]*M[4,3]-M[1,1]*M[2,4]*M[3,3]*M[4,2]-M[1,2]*M[2,1]*M[3,3]*M[4,4]+M[1,2]*M[2,1]*M[3,4]*M[4,3]+M[1,2]*M[2,3]*M[3,1]*M[4,4]-M[1,2]*M[2,3]*M[3,4]*M[4,1]-M[1,2]*M[2,4]*M[3,1]*M[4,3]+M[1,2]*M[2,4]*M[3,3]*M[4,1]+M[1,3]*M[2,1]*M[3,2]*M[4,4]-M[1,3]*M[2,1]*M[3,4]*M[4,2]-M[1,3]*M[2,2]*M[3,1]*M[4,4]+M[1,3]*M[2,2]*M[3,4]*M[4,1]+M[1,3]*M[2,4]*M[3,1]*M[4,2]-M[1,3]*M[2,4]*M[3,2]*M[4,1]-M[1,4]*M[2,1]*M[3,2]*M[4,3]+M[1,4]*M[2,1]*M[3,3]*M[4,2]+M[1,4]*M[2,2]*M[3,1]*M[4,3]-M[1,4]*M[2,2]*M[3,3]*M[4,1]-M[1,4]*M[2,3]*M[3,1]*M[4,2]+M[1,4]*M[2,3]*M[3,2]*M[4,1]
end

function D1_noH(M)
  return M[1, 1] + M[1, 2] + M[2, 1] + M[2, 2]
end

function D2_noH(M)
  return M[1,1]^2+2*M[1,1]*M[2,1]-M[1,2]^2-2*M[1,2]*M[2,2]-M[1,3]^2-2*M[1,3]*M[2,3]-M[1,4]^2-2*M[1,4]*M[2,4]+M[2,1]^2-M[2,2]^2-M[2,3]^2-M[2,4]^2
end

function D3_noH(M)
  return M[1,1]^3+(M[1,2]+M[2,1]-M[2,2])*M[1,1]^2+(-M[1,2]^2+(2*M[2,1]-2*M[2,2])*M[1,2]-M[3,2]^2-2*M[3,1]*M[3,2]-M[3,3]^2+2*M[3,3]*M[4,4]-M[3,4]^2-2*M[3,4]*M[4,3]-M[4,1]^2-2*M[4,1]*M[4,2]-M[4,2]^2-M[4,3]^2-M[4,4]^2-M[1,3]^2-2*M[1,3]*M[2,3]-M[1,4]^2-2*M[1,4]*M[2,4]-M[2,1]^2-2*M[2,1]*M[2,2]-M[2,2]^2-M[2,3]^2-M[2,4]^2-M[3,1]^2)*M[1,1]-M[1,2]^3+(M[2,1]-M[2,2])*M[1,2]^2+(-M[1,3]^2-2*M[1,3]*M[2,3]-M[1,4]^2-2*M[1,4]*M[2,4]+M[2,1]^2+2*M[2,1]*M[2,2]+M[2,2]^2-M[2,3]^2-M[2,4]^2+M[3,1]^2+2*M[3,1]*M[3,2]+M[3,2]^2-M[3,3]^2+2*M[3,3]*M[4,4]-M[3,4]^2-2*M[3,4]*M[4,3]+M[4,1]^2+2*M[4,1]*M[4,2]+M[4,2]^2-M[4,3]^2-M[4,4]^2)*M[1,2]-M[2,1]^3-M[2,1]^2*M[2,2]+(M[1,3]^2+2*M[1,3]*M[2,3]+M[1,4]^2+2*M[1,4]*M[2,4]+M[2,2]^2+M[2,3]^2+M[2,4]^2-M[3,1]^2-2*M[3,1]*M[3,2]-M[3,2]^2-M[3,3]^2+2*M[3,3]*M[4,4]-M[3,4]^2-2*M[3,4]*M[4,3]-M[4,1]^2-2*M[4,1]*M[4,2]-M[4,2]^2-M[4,3]^2-M[4,4]^2)*M[2,1]+M[2,2]^3+(M[1,3]^2+2*M[1,3]*M[2,3]+M[1,4]^2+2*M[1,4]*M[2,4]+M[2,3]^2+M[2,4]^2+M[3,1]^2+2*M[3,1]*M[3,2]+M[3,2]^2-M[3,3]^2+2*M[3,3]*M[4,4]-M[3,4]^2-2*M[3,4]*M[4,3]+M[4,1]^2+2*M[4,1]*M[4,2]+M[4,2]^2-M[4,3]^2-M[4,4]^2)*M[2,2]+((2*M[3,3]-2*M[4,4])*M[3,1]+(2*M[3,3]-2*M[4,4])*M[3,2]+2*(M[4,1]+M[4,2])*(M[3,4]+M[4,3]))*M[1,3]+((2*M[3,4]+2*M[4,3])*M[3,1]+(2*M[3,4]+2*M[4,3])*M[3,2]-2*(M[4,1]+M[4,2])*(M[3,3]-M[4,4]))*M[1,4]+((2*M[3,3]-2*M[4,4])*M[3,1]+(2*M[3,3]-2*M[4,4])*M[3,2]+2*(M[4,1]+M[4,2])*(M[3,4]+M[4,3]))*M[2,3]+2*M[2,4]*((M[3,4]+M[4,3])*M[3,1]+(M[3,4]+M[4,3])*M[3,2]-(M[4,1]+M[4,2])*(M[3,3]-M[4,4]))
end

function D4_noH(M)
  return M[1,1]^4+M[1,2]^4+M[1,3]^4+M[2,1]^4+M[2,2]^4+M[2,3]^4+M[3,1]^4+M[3,2]^4+M[3,3]^4+M[1,4]^4+M[2,4]^4+M[3,4]^4+8*M[3,3]*M[3,4]*M[4,3]*M[4,4]-8*M[4,1]*(M[3,2]*M[4,2]+M[3,3]*M[4,3]+M[3,4]*M[4,4])*M[3,1]+8*M[4,2]*(M[3,3]*M[4,3]+M[3,4]*M[4,4])*M[3,2]+8*M[2,4]*(M[3,3]*M[3,4]+M[4,3]*M[4,4])*M[2,3]+(M[4,1]^2-M[4,2]^2-M[4,3]^2-M[4,4]^2)^2+(2*M[4,1]^2-2*M[4,2]^2-2*M[4,3]^2+2*M[4,4]^2)*M[3,4]^2+(2*M[3,4]^2+2*M[4,1]^2-2*M[4,2]^2+2*M[4,3]^2-2*M[4,4]^2)*M[3,3]^2+(2*M[3,3]^2+2*M[3,4]^2+2*M[4,1]^2+2*M[4,2]^2-2*M[4,3]^2-2*M[4,4]^2)*M[3,2]^2+(-2*M[3,2]^2-2*M[3,3]^2-2*M[3,4]^2+2*M[4,1]^2+2*M[4,2]^2+2*M[4,3]^2+2*M[4,4]^2)*M[3,1]^2+(2*M[3,1]^2-2*M[3,2]^2-2*M[3,3]^2+2*M[3,4]^2+2*M[4,1]^2-2*M[4,2]^2-2*M[4,3]^2+2*M[4,4]^2)*M[2,4]^2+(2*M[2,4]^2+2*M[3,1]^2-2*M[3,2]^2+2*M[3,3]^2-2*M[3,4]^2+2*M[4,1]^2-2*M[4,2]^2+2*M[4,3]^2-2*M[4,4]^2)*M[2,3]^2+((8*M[3,2]*M[3,3]+8*M[4,2]*M[4,3])*M[2,3]+8*M[2,4]*(M[3,2]*M[3,4]+M[4,2]*M[4,4]))*M[2,2]+(2*M[2,3]^2+2*M[2,4]^2+2*M[3,1]^2+2*M[3,2]^2-2*M[3,3]^2-2*M[3,4]^2+2*M[4,1]^2+2*M[4,2]^2-2*M[4,3]^2-2*M[4,4]^2)*M[2,2]^2+((-8*M[3,1]*M[3,2]-8*M[4,1]*M[4,2])*M[2,2]+(-8*M[3,1]*M[3,3]-8*M[4,1]*M[4,3])*M[2,3]-8*M[2,4]*(M[3,1]*M[3,4]+M[4,1]*M[4,4]))*M[2,1]+(-2*M[2,2]^2-2*M[2,3]^2-2*M[2,4]^2+2*M[3,1]^2+2*M[3,2]^2+2*M[3,3]^2+2*M[3,4]^2+2*M[4,1]^2+2*M[4,2]^2+2*M[4,3]^2+2*M[4,4]^2)*M[2,1]^2+((-8*M[3,2]*M[4,3]+8*M[3,3]*M[4,2])*M[2,1]+(8*M[3,1]*M[4,3]-8*M[3,3]*M[4,1])*M[2,2]+8*M[2,3]*(-M[3,1]*M[4,2]+M[3,2]*M[4,1]))*M[1,4]+(-2*M[2,1]^2+2*M[2,2]^2+2*M[2,3]^2-2*M[2,4]^2-2*M[3,1]^2+2*M[3,2]^2+2*M[3,3]^2-2*M[3,4]^2-2*M[4,1]^2+2*M[4,2]^2+2*M[4,3]^2-2*M[4,4]^2)*M[1,4]^2+((-8*M[2,3]*M[2,4]-8*M[3,3]*M[3,4]-8*M[4,3]*M[4,4])*M[1,4]+(8*M[3,2]*M[4,4]-8*M[3,4]*M[4,2])*M[2,1]+(-8*M[3,1]*M[4,4]+8*M[3,4]*M[4,1])*M[2,2]-8*M[2,4]*(-M[3,1]*M[4,2]+M[3,2]*M[4,1]))*M[1,3]+(2*M[1,4]^2-2*M[2,1]^2+2*M[2,2]^2-2*M[2,3]^2+2*M[2,4]^2-2*M[3,1]^2+2*M[3,2]^2-2*M[3,3]^2+2*M[3,4]^2-2*M[4,1]^2+2*M[4,2]^2-2*M[4,3]^2+2*M[4,4]^2)*M[1,3]^2+((-8*M[2,2]*M[2,3]-8*M[3,2]*M[3,3]-8*M[4,2]*M[4,3])*M[1,3]+(-8*M[2,2]*M[2,4]-8*M[3,2]*M[3,4]-8*M[4,2]*M[4,4])*M[1,4]+(-8*M[3,3]*M[4,4]+8*M[3,4]*M[4,3])*M[2,1]+(8*M[3,1]*M[4,4]-8*M[3,4]*M[4,1])*M[2,3]+8*M[2,4]*(-M[3,1]*M[4,3]+M[3,3]*M[4,1]))*M[1,2]+(2*M[1,3]^2+2*M[1,4]^2-2*M[2,1]^2-2*M[2,2]^2+2*M[2,3]^2+2*M[2,4]^2-2*M[3,1]^2-2*M[3,2]^2+2*M[3,3]^2+2*M[3,4]^2-2*M[4,1]^2-2*M[4,2]^2+2*M[4,3]^2+2*M[4,4]^2)*M[1,2]^2+((8*M[2,1]*M[2,2]+8*M[3,1]*M[3,2]+8*M[4,1]*M[4,2])*M[1,2]+(8*M[2,1]*M[2,3]+8*M[3,1]*M[3,3]+8*M[4,1]*M[4,3])*M[1,3]+(8*M[2,1]*M[2,4]+8*M[3,1]*M[3,4]+8*M[4,1]*M[4,4])*M[1,4]+(8*M[3,3]*M[4,4]-8*M[3,4]*M[4,3])*M[2,2]+(-8*M[3,2]*M[4,4]+8*M[3,4]*M[4,2])*M[2,3]+8*M[2,4]*(M[3,2]*M[4,3]-M[3,3]*M[4,2]))*M[1,1]+(-2*M[1,2]^2-2*M[1,3]^2-2*M[1,4]^2-2*M[2,1]^2-2*M[2,2]^2-2*M[2,3]^2-2*M[2,4]^2-2*M[3,1]^2-2*M[3,2]^2-2*M[3,3]^2-2*M[3,4]^2-2*M[4,1]^2-2*M[4,2]^2-2*M[4,3]^2-2*M[4,4]^2)*M[1,1]^2
end

# Sylvester Criterion
function SC(M)
  d1 = M[1,1]
  d2 = det(M[1:2,1:2])
  d3 = det(M[1:3,1:3])
  d4 = det(M)
  return ((real(d1) > 0) & (real(d2) > 0) & (real(d3) > 0) & (real(d4) > 0))  
end

function SCdec(M)
  return ((real(D1(M)) > 0) & (real(D2(M)) > 0) & (real(D3(M)) > 0) & (real(D4(M)) > 0))
end

function SCdec_noH(M)
  return ((D1_noH(M) > 0) & (D2_noH(M) > 0) & (D3_noH(M) > 0) & (D4_noH(M) > 0))
end

println("Time to apply the function H to the array of matrices:")
# Applying the function H to the array of matrices
@btime global HM = [H(MMlist[s]) for s in 1:sz];
println()


println("Time to check positive definite using Julia's built-in function (Choletsky), parallel:")
# Checking positive definite using characteristic polynomial
num_threads = Threads.nthreads()
sz_threaded = Int32(sz / num_threads)
@btime global HMosdef_noH_par = Folds.collect([isposdef(HM[s]) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()


println("Time to check positive definite using Julia's built-in function (Choletsky):")
# Checking positive definite using Julia's built-in function (Choletsky)
@btime global HMposdef = [isposdef(HM[s]) for s in 1:sz];
println()

println("Time to check positive definite using Julia's built-in function (Choletsky), avoiding H, parallel:")
# Checking positive definite using characteristic polynomial
num_threads = Threads.nthreads()
sz_threaded = Int32(sz / num_threads)
@btime global HMosdef_noH_par = Folds.collect([isposdef(H(MMlist[s])) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()

println("Time to check positive definite using Julia's built-in function (Choletsky), avoiding H:")
# Checking positive definite using Julia's built-in function (Choletsky)
@btime global HMposdef = [isposdef(H(MMlist[s])) for s in 1:sz];
println()


# Creating empty arrays to store the non-physical matrices
nonposdef = Array{Int32}(undef, sz)
nonpossemdef = Array{Int32}(undef, sz)
nonposdefSC = Array{Int32}(undef, sz)
nonpossemdefCP = Array{Int32}(undef, sz)
nonposdefSC_noH = Array{Int32}(undef, sz)
nonpossemdefCP_noH = Array{Int32}(undef, sz)
nonpossemdefCP_noH_par = Array{Int32}(undef, sz)

println("Time to check positive definite using Sylvester's Criterion:")
# Checking positive definite using Sylvester's Criterion
@btime global SCHMposdef = [SC(HM[s]) for s in 1:sz];
println()


println("Time to check positive definite using Sylvester's Criterion with determinants calculated as separate functions in Julia:")
# Checking positive definite using Sylvester's Criterion
@btime global SCHMposdef = [SCdec(HM[s]) for s in 1:sz];
println()


println("Time to check positive definite using Sylvester's Criterion with determinants calculated as separate functions in Julia avoiding H, parallel:")
# Checking positive definite using characteristic polynomial
num_threads = Threads.nthreads()
sz_threaded = Int32(sz / num_threads)
@btime global SCHMpossemdef_noH_par = Folds.collect([SCdec_noH(MMlist[s]) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()

println("Time to check positive definite using Sylvester's Criterion with determinants calculated as separate functions in Julia avoiding H:")
# Checking positive definite using Sylvester's Criterion
@btime global SCHMposdef_noH = [SCdec_noH(MMlist[s]) for s in 1:sz];
println()

println("Time to check positive semi-definite using characteristic polynomial:")
# Checking positive definite using characteristic polynomial
@btime global CPHMpossemdef = [PSD(HM[s]) for s in 1:sz];
println()


println("Time to check positive semi-definite using characteristic polynomial avoiding H, parallel:")
# Checking positive definite using characteristic polynomial
num_threads = Threads.nthreads()
sz_threaded = Int32(sz / num_threads)
@btime global CPHMpossemdef_noH_par = Folds.collect([PSD_noH(MMlist[s]) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()

println("Time to check positive semi-definite using characteristic polynomial avoiding H:")
# Checking positive definite using characteristic polynomial
@btime global CPHMpossemdef_noH = [PSD_noH(MMlist[s]) for s in 1:sz];
println()

println("Time to check positive semi-definite using characteristic polynomial avoiding H in a for-loop:")
# Checking positive definite using characteristic polynomial
@btime @inbounds for s in 1:sz @inbounds CPHMpossemdef_noH[s] = PSD_noH(MMlist[s]) end 
println()

println("Time to check positive semi-definite using characteristic polynomial avoiding H in a for-loop, parallel:")
# Checking positive definite using characteristic polynomial
@btime Threads.@threads for s in 1:num_threads @inbounds for ss in ((s - 1) * sz_threaded + 1):(s * sz_threaded) @inbounds CPHMpossemdef_noH[ss] = PSD_noH(MMlist[ss]) end end
println()

println("Time to check positive definite in a for-loop using built-in function:")
# Checking positive definite in a for-loop using built-in function
i0 = 0
@time for s in 1:sz
  global i0
  if isposdef(HM[s]) == false 
    i0 += 1
    nonposdef[i0] = s
  end
end
println()

println("The number of non-positive definite matrices (according to Choletsky) is ", i0, ".")
println()
println("The percentage of non-positive definite matrices (according to Choletsky) is ", round(i0/sz*100, digits = 3), "%.")
println()

i1 = 0
for s in 1:sz
  global i1
  if SCHMposdef[s] == false
    i1 += 1
    nonposdefSC[i1] = s
  end
end
println()

println("The number of non-positive definite matrices (according to Sylvester) is ", i1, ".")
println()

i11 = 0
for s in 1:sz
  global i11
  if SCHMposdef_noH[s] == false
    i11 += 1
    nonposdefSC_noH[i11] = s
  end
end
println()

println("The number of non-positive definite matrices (according to Sylvester, avoiding H) is ", i11, ".")
println()

i2 = 0
for s in 1:sz
  global i2
  if CPHMpossemdef[s] == false
    i2 += 1
    nonpossemdefCP[i2] = s
  end
end

println("The number of non-positive semi-definite matrices (according to characteristic polynomial) is ", i2, ".")
println()

i21 = 0
for s in 1:sz
  global i21
  if CPHMpossemdef_noH[s] == false
    i21 += 1
    nonpossemdefCP_noH[i21] = s
  end
end

println("The number of non-positive semi-definite matrices (according to characteristic polynomial, avoiding H) is ", i21, ".")
println()

i22 = 0
for s in 1:num_threads
  for j in 1:sz_threaded 
    global i22
    if CPHMpossemdef_noH_par[s][j] == false
      i22 += 1
      nonpossemdefCP_noH_par[i22] = s
    end
  end
end

println("The number of non-positive semi-definite matrices (according to characteristic polynomial, avoiding H, parallel) is ", i22, ".")
println()


println("The time to check positive semi-definite by computing eigenvalues, parallel:")
num_threads = Threads.nthreads()
sz_threaded = Int32(sz / num_threads)
@btime global HMev_par = Folds.collect([all(eigvals(HM[s]) .>= 0) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()

println("The time to check positive semi-definite by computing eigenvalues:")
# Checking positive semi-definite by computing eigenvalues
@btime global HMev = [all(eigvals(HM[s]) .>= 0) for s in 1:sz];
println()

println("The time to check positive semi-definite by computing eigenvalues, avoiding H, parallel:")
num_threads = Threads.nthreads()
sz_threaded = Int32(sz / num_threads)
@btime global HMev_par = Folds.collect([all(eigvals(H(MMlist[s])) .>= 0) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)]  for i in 1:num_threads);
println()

println("The time to check positive semi-definite by computing eigenvalues, avoiding H:")
# Checking positive semi-definite by computing eigenvalues
@btime global HMev = [all(eigvals(H(MMlist[s])) .>= 0) for s in 1:sz];
println()


println("The time to create an array of all eigenvalues for all matrices:")
@btime global HMev = [eigvals(HM[s]) for s in 1:sz];
println()

# Counting the number of not positive-semidefinite matrices using eigenvalues
evnegnum = 0
for s in 1:sz
  if eigvals(HM[s])[1] < 0
    global evnegnum += 1
    nonpossemdef[evnegnum] = s
  end 
end

println("The number of matrices that are not positive semi-definite matrices (according to eigenvalue computation) is ", evnegnum)
println()
println("Checking if the indices of non positive definite (by Choletsky) are the same as the indices of not positive semi-definite (according to eigenvalues):")
println(nonposdef[1:i0] == nonpossemdef[1:evnegnum])
println()
println("Checking if the indices of non positive definite (by Choletsky) are the same as the indices of non positive definite (Sylvester Criterion):")
println(nonposdef[1:i0] == nonposdefSC[1:i1])
println()
println("Checking if the indices of non positive definite (by Choletsky) are the same as the indices of non positive semi-definite (characteristic polynomial):")
println(nonposdef[1:i0] == nonpossemdefCP[1:i1])
println()
