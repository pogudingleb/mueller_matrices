using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using Folds

# loading external files with matrix computation
include("PSD.jl")
include("SC.jl")
include("H.jl")

# Size
sz = 600*700

# Setup for multi-threading
num_threads = Threads.nthreads()
sz_threaded = Int32(sz / num_threads)

# Download the input data into a matrix of arrays
println()
println("=================================")
println()
println("Benchmarking the performance of verifying positive definite and positive semidefinite properties of Mueller matrices using Choletsky, Sylvester's criterion, characteristic polynomial, and eigenvalues. This is done using @btime, which can take some time to run.")
println()

println("=================================")
println()
println("The CPU in this system is ", Sys.cpu_info()[1].model, ".")
println()
println("The total number of CPU cores is ", Sys.CPU_THREADS, ".")
println()
println("The number of used threads in this Julia run is ", num_threads, ".")
println()


println("=================================")
println()
println("Time to:")
println()
println("- download the input data into a matrix of arrays (the way the data are currently stored):")
@btime global MM = [reshape(readdlm("../data/MM_"*string(i)*string(j)*".dat"), 1, sz) for i in 1:4, j in 1:4];
println()

# Convert the matrix of arrays into an array of matrices

println("- convert the matrix of arrays into an array of matrices:")
@btime global MMlist = [[MM[i,j][s] for i in 1:4, j in 1:4] for s in 1:sz];
println()


println("- apply the function H to the array of matrices:")
# Applying the function H to the array of uploaded matrices
@btime global HM = [H(MMlist[s]) for s in 1:sz];
println()


println("===================================")
println()
println("Time to check positive:")
println()
println("- definite using Julia's built-in function (Cholesky), parallel:")
@btime global HMosdef_noH_par = Folds.collect([isposdef(HM[s]) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()


println("- definite using Julia's built-in function (Cholesky):")
@btime global HMposdef = [isposdef(HM[s]) for s in 1:sz];
println()

println("- definite using Julia's built-in function (Cholesky), not precomputing H, parallel:")
@btime global HMosdef_noH_par = Folds.collect([isposdef(H(MMlist[s])) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()

println("- definite using Julia's built-in function (Cholesky), not precomputing H:")
@btime global HMposdef_noH = [isposdef(H(MMlist[s])) for s in 1:sz];
println()


println("- definite using Sylvester's Criterion, not precomputing H, parallel:")
@btime global SCHMpossemdef_noH_par = Folds.collect([SC(MMlist[s]) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()

println("- definite using Sylvester's Criterion, not precomputing H:")
@btime global SCHMposdef_noH = [SC(MMlist[s]) for s in 1:sz];
println()


println("- semi-definite using characteristic polynomial, not precomputing H, parallel:")
@btime global CPHMpossemdef_noH_par = Folds.collect([PSD(MMlist[s]) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()

println("- semi-definite using characteristic polynomial, not precomputing H:")
@btime global CPHMpossemdef_noH = [PSD(MMlist[s]) for s in 1:sz];
println()


println("- semi-definite by computing eigenvalues, parallel:")
@btime global HMev_par = Folds.collect([all(eigvals(HM[s]) .>= 0) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)] for i in 1:num_threads);
println()

println("- semi-definite by computing eigenvalues:")
@btime global HMev = [all(eigvals(HM[s]) .>= 0) for s in 1:sz];
println()

println("- semi-definite by computing eigenvalues, not precomputing H, parallel:")
@btime global HMev_par = Folds.collect([all(eigvals(H(MMlist[s])) .>= 0) for s in ((i - 1) * sz_threaded + 1):(i * sz_threaded)]  for i in 1:num_threads);
println()

println("- semi-definite by computing eigenvalues, not precomputing H:")
@btime global HMev = [all(eigvals(H(MMlist[s])) .>= 0) for s in 1:sz];
println()

println("=================================")
println()
println("The number of non-positive:")
println()

# Storing indices of non-positive definite matrices in a for-loop using built-in Cholesky
Cholesky_count = 0
nonposdef = Array{Int32}(undef, sz)
for s in 1:sz
  global Cholesky_count
  if isposdef(HM[s]) == false 
    Cholesky_count += 1
    nonposdef[Cholesky_count] = s
  end
end

println("- definite matrices (according to Choletsky) is ", Cholesky_count)
println("-     (the percentage of non-positive definite matrices is ", round(Cholesky_count/sz*100, digits = 3), "%.")
println()

nonposdefSC = Array{Int32}(undef, sz)
Sylvester_count = 0
for s in 1:sz
  global Sylvester_count
  if SCHMposdef_noH[s] == false
    Sylvester_count += 1
    nonposdefSC[Sylvester_count] = s
  end
end

println("- definite matrices (according to Sylvester) is ", Sylvester_count, ".")
println()

nonpossemdefCP = Array{Int32}(undef, sz)
CharPoly_count = 0
for s in 1:sz
  global CharPoly_count
  if CPHMpossemdef_noH[s] == false
    CharPoly_count += 1
    nonpossemdefCP[CharPoly_count] = s
  end
end

println("- semi-definite matrices (according to characteristic polynomial) is ", CharPoly_count, ".")
println()

# Counting the number of not positive-semidefinite matrices using eigenvalues
nonpossemdef = Array{Int32}(undef, sz)
evnegnum = 0
for s in 1:sz
  if eigvals(HM[s])[1] < 0
    global evnegnum += 1
    nonpossemdef[evnegnum] = s
  end 
end

println("- semi-definite matrices (according to eigenvalue computation) is ", evnegnum)
println()
println("=================================")
println()
println("Checking if the indices of non positive definite (by Cholesky) are the same as the indices of not positive")
println()
println("- semi-definite (according to eigenvalues): ",
         nonposdef[1:Cholesky_count] == nonpossemdef[1:evnegnum])
println()
println("- definite (Sylvester Criterion): ", 
        nonposdef[1:Cholesky_count] == nonposdefSC[1:Sylvester_count])
println()
println("- semi-definite (characteristic polynomial): ", 
        nonposdef[1:Cholesky_count] == nonpossemdefCP[1:CharPoly_count])
