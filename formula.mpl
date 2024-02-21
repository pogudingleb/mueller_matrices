with(LinearAlgebra):

# This H is in fact 1/2 * H in terms of the paper, the scaling is taken to reduce the numer of constants in the expressions
# Therefore, the formulas in the paper differ by scaling by an appropriate power of 1/2
H := 1/2 * Matrix([
    [M[1,1] + M[1,2] + M[2,1] + M[2,2], M[1,3] + M[2,3] + I * (M[1,4] + M[2,4]), M[3,1] + M[3,2] - I * (M[4,1] + M[4,2]), M[3,3] + M[4,4] + I * (M[3,4] - M[4,3])],
    [M[1,3] + M[2,3] - I * (M[1,4] + M[2,4]), M[1,1] - M[1,2] + M[2,1] - M[2,2], M[3,3] - M[4,4] - I * (M[3,4] + M[4,3]), M[3,1] - M[3,2] - I * (M[4,1] - M[4,2])],
    [M[3,1] + M[3,2] + I * (M[4,1] + M[4,2]), M[3,3] - M[4,4] + I * (M[3,4] + M[4,3]), M[1,1] + M[1,2] - M[2,1] - M[2,2], M[1,3] - M[2,3] + I * (M[1,4] - M[2,4])],
    [M[3,3] + M[4,4] - I * (M[3,4] - M[4,3]), M[3,1] - M[3,2] + I * (M[4,1] - M[4,2]), M[1,3] - M[2,3] - I * (M[1,4] - M[2,4]), M[1,1] - M[1,2] - M[2,1] + M[2,2]]
]);

tr1 := Trace(H);
tr2 := expand(Trace(H^2));
tr3 := expand(Trace(H^3));
tr4 := expand(Trace(H^4));

#####################
# The formula for tr3
#####################

prod1x2xs := M[1, 2] * M[2, 2] + M[1, 3] * M[2, 3] + M[1, 4] * M[2, 4];
prod1x3xs := M[1, 2] * M[3, 2] + M[1, 3] * M[3, 3] + M[1, 4] * M[3, 4];
prod1x4xs := M[1, 2] * M[4, 2] + M[1, 3] * M[4, 3] + M[1, 4] * M[4, 4];
with_pairs_mixed := M[2, 1] * prod1x2xs + M[3, 1] * prod1x3xs + M[4, 1] * prod1x4xs;
corner_det := Determinant([[M[2, 2], M[2, 3], M[2, 4]], [M[3, 2], M[3, 3], M[3, 4]], [M[4, 2], M[4, 3], M[4, 4]]]); # det(\tilde{M}) in the notation of the paper

print("Verification of the formula for tr(H^3): ", expand(tr3 - 3/2 * M[1, 1] * tr2 + M[1, 1]^3 - 3 * (with_pairs_mixed + corner_det)));

#####################
# The formula for tr3
#####################

# det(M) in the paper
full_det := Determinant(Matrix([seq([seq(M[i, j], j=1..4)], i=1..4)]));

# S_1, S_2, S_3 in the paper
row_squaresums := [seq(add([seq(M[i, j]^2, j=2..4)]), i=2..4)];

# A in the paper
row_squaresums_sq := expand(add([seq(row_squaresums[i]^2, i=1..3)]));

# B in the paper
anticommuting_with1 := expand(M[2, 1]^2 * (row_squaresums[2] + row_squaresums[3]) + M[3, 1]^2 * (row_squaresums[1] + row_squaresums[3]) + M[4, 1]^2 * (row_squaresums[1] + row_squaresums[2]));

# P_{i - 1, j - 1} in the paper
prod1x2x := M[1, 1] * M[2, 1] + prod1x2xs;
prod1x3x := M[1, 1] * M[3, 1] + prod1x3xs;
prod1x4x := M[1, 1] * M[4, 1] + prod1x4xs;
prod2x3x := M[2, 1] * M[3, 1] - M[2, 2] * M[3, 2] - M[2, 3] * M[3, 3] - M[2, 4] * M[3, 4];
prod2x4x := M[2, 1] * M[4, 1] - M[2, 2] * M[4, 2] - M[2, 3] * M[4, 3] - M[2, 4] * M[4, 4];
prod3x4x := M[3, 1] * M[4, 1] - M[3, 2] * M[4, 2] - M[3, 3] * M[4, 3] - M[3, 4] * M[4, 4];

# C in the paper
halfmixed_products := expand(-prod1x2x^2 - prod1x3x^2 - prod1x4x^2 + prod2x3x^2 + prod2x4x^2 + prod3x4x^2);

# D in the paper
squares1x := add([seq(M[1, i]^2, i=2..4)]);

# F in the paper
squaresx1 := add([seq(M[i, 1]^2, i=2..4)]);

# S_1 + S_2 + S_3
corner_squaresum := add(row_squaresums);

print("Verification of the formula for tr(H^4): ", expand(tr4 - (-2 * full_det + 3/4 * (tr2)^2 - 1/2 * row_squaresums_sq - anticommuting_with1 - squares1x * (1/2 * squares1x + corner_squaresum) + 4 * M[1, 1] * (2 * corner_det + with_pairs_mixed) - halfmixed_products - M[1, 1]^2 * squaresx1 - 1/2 * (M[1, 1]^4 + M[2, 1]^4 + M[3, 1]^4 + M[4, 1]^4))));
