with(LinearAlgebra):

H := 1/2 * Matrix([
    [M[1,1] + M[1,2] + M[2,1] + M[2,2], M[1,3] + M[2,3] + I * (M[1,4] + M[2,4]), M[3,1] + M[3,2] - I * (M[4,1] + M[4,2]), M[3,3] + M[4,4] + I * (M[3,4] - M[4,3])],
    [M[1,3] + M[2,3] - I * (M[1,4] + M[2,4]), M[1,1] - M[1,2] + M[2,1] - M[2,2], M[3,3] - M[4,4] - I * (M[3,4] + M[4,3]), M[3,1] - M[3,2] - I * (M[4,1] - M[4,2])],
    [M[3,1] + M[3,2] + I * (M[4,1] + M[4,2]), M[3,3] - M[4,4] + I * (M[3,4] + M[4,3]), M[1,1] + M[1,2] - M[2,1] - M[2,2], M[1,3] - M[2,3] + I * (M[1,4] - M[2,4])],
    [M[3,3] + M[4,4] - I * (M[3,4] - M[4,3]), M[3,1] - M[3,2] + I * (M[4,1] - M[4,2]), M[1,3] - M[2,3] - I * (M[1,4] - M[2,4]), M[1,1] - M[1,2] - M[2,1] + M[2,2]]
]);

tr1 := Trace(H);
tr2 := expand(Trace(H^2)); # 16 mult
tr3 := expand(Trace(H^3));
tr4 := expand(Trace(H^4));

# The formula for tr3

prod1x2xs := M[1, 2] * M[2, 2] + M[1, 3] * M[2, 3] + M[1, 4] * M[2, 4]; # 3 mult
prod1x3xs := M[1, 2] * M[3, 2] + M[1, 3] * M[3, 3] + M[1, 4] * M[3, 4]; # 3 mult
prod1x4xs := M[1, 2] * M[4, 2] + M[1, 3] * M[4, 3] + M[1, 4] * M[4, 4]; # 3 mult
with_pairs_mixed := M[2, 1] * prod1x2xs + M[3, 1] * prod1x3xs + M[4, 1] * prod1x4xs; # 3 mult
corner_det := Determinant([[M[2, 2], M[2, 3], M[2, 4]], [M[3, 2], M[3, 3], M[3, 4]], [M[4, 2], M[4, 3], M[4, 4]]]); # 9 mult

print("Verification of the formula for tr(H^3): ", expand(tr3 - 3/2 * M[1, 1] * tr2 + M[1, 1]^3 - 3 * (with_pairs_mixed + corner_det))); # 2 mult

# for tr3 --> 23 mult

full_det := Determinant(Matrix([seq([seq(M[i, j], j=1..4)], i=1..4)]));

tr4_1 := expand(tr4 + 2 * full_det);

tr4_2 := expand(tr4_1 - 3 / 4 * (tr2)^2);

row_squaresums := [seq(add([seq(M[i, j]^2, j=2..4)]), i=2..4)];
row_squaresums_sq := expand(add([seq(row_squaresums[i]^2, i=1..3)]));

tr4_3 := expand(tr4_2 + 1/2 * row_squaresums_sq);

anticommuting_with1 := expand(M[2, 1]^2 * (row_squaresums[2] + row_squaresums[3]) + M[3, 1]^2 * (row_squaresums[1] + row_squaresums[3]) + M[4, 1]^2 * (row_squaresums[1] + row_squaresums[2]));

tr4_4 := tr4_3 + anticommuting_with1;

squares1x := add([seq(M[1, i]^2, i=2..4)]);
squaresx1 := add([seq(M[i, 1]^2, i=2..4)]);

tr4_5 := tr4_4 + 1/2 * expand(squares1x^2 + squaresx1^2);

tr4_6 := tr4_5 - 8 * expand(M[1, 1] * corner_det);

tr4_7 := tr4_6 - 4 * expand(M[1, 1] * with_pairs_mixed);

prod1x2x := M[1, 1] * M[2, 1] + prod1x2xs;
prod1x3x := M[1, 1] * M[3, 1] + prod1x3xs;
prod1x4x := M[1, 1] * M[4, 1] + prod1x4xs;
prod2x3x := M[2, 1] * M[3, 1] - M[2, 2] * M[3, 2] - M[2, 3] * M[3, 3] - M[2, 4] * M[3, 4];
prod2x4x := M[2, 1] * M[4, 1] - M[2, 2] * M[4, 2] - M[2, 3] * M[4, 3] - M[2, 4] * M[4, 4];
prod3x4x := M[3, 1] * M[4, 1] - M[3, 2] * M[4, 2] - M[3, 3] * M[4, 3] - M[3, 4] * M[4, 4];
halfmixed_products := expand(-prod1x2x^2 - prod1x3x^2 - prod1x4x^2 + prod2x3x^2 + prod2x4x^2 + prod3x4x^2);
tr4_8 := expand(tr4_7 + halfmixed_products);

corner_squaresum := add(row_squaresums);
products_two_ones := expand(corner_squaresum * squares1x);
products_four_ones :=expand((M[1, 1]^2 - squaresx1)^2);
tr4_9 := tr4_8 + products_two_ones - 1/2 * products_four_ones;
tr4_10 := tr4_9 + M[1, 1]^4 + 1/2 * (M[2, 1]^4 + M[3, 1]^4 + M[4, 1]^4);
