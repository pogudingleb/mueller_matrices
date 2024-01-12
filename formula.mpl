with(LinearAlgebra):

H := Matrix([
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

F2 := M[2, 2] * M[1, 2] + M[2, 3] * M[1, 3] + M[2, 4] * M[1, 4]; # 3 mult
F3 := M[3, 2] * M[1, 2] + M[3, 3] * M[1, 3] + M[3, 4] * M[1, 4]; # 3 mult
F4 := M[4, 2] * M[1, 2] + M[4, 3] * M[1, 3] + M[4, 4] * M[1, 4]; # 3 mult
with_pairs_mixed := M[2, 1] * F2 + M[3, 1] * F3 + M[4, 1] * F4; # 3 mult
corner_det := Determinant([[M[2, 2], M[2, 3], M[2, 4]], [M[3, 2], M[3, 3], M[3, 4]], [M[4, 2], M[4, 3], M[4, 4]]]); # 9 mult

print("Verification of the formula for tr(H^3): ", expand(tr3 - 3 * M[1, 1] * tr2 + 8 * M[1, 1]^3 - 24 * (with_pairs_mixed + corner_det))); # 2 mult

# for tr3 --> 23 mult

full_det := Determinant(Matrix([seq([seq(M[i, j], j=1..4)], i=1..4)]));

tr4_1 := expand(tr4 + 32 * full_det);

tr4_2 := expand(tr4_1 - 12 * (tr2 / 4)^2);

G_tail := [seq(add([seq(M[i, j]^2, j=2..4)]), i=2..4)];
G_tail2 := expand(add([seq(G_tail[i]^2, i=1..3)]));

tr4_3 := expand(tr4_2 + 8 * G_tail2);

tr4_4 := tr4_3;

anticommuting_with1_op := expand(M[2, 1]^2 * (G_tail[2] + G_tail[3]) + M[3, 1]^2 * (G_tail[1] + G_tail[3]) + M[4, 1]^2 * (G_tail[1] + G_tail[2]));

tr4_5 := tr4_4 + 16 * anticommuting_with1_op;

squares1x := add([seq(M[1, i]^2, i=2..4)]);
squaresx1 := add([seq(M[i, 1]^2, i=2..4)]);

tr4_6 := tr4_5 + 8 * expand(squares1x^2 + squaresx1^2);

tr4_7 := tr4_6 - 128 * expand(M[1, 1] * corner_det);

tr4_8 := tr4_7 - 64 * expand(M[1, 1] * with_pairs_mixed);

S12 := M[1, 1] * M[2, 1] + M[1, 2] * M[2, 2] + M[1, 3] * M[2, 3] + M[1, 4] * M[2, 4];
S13 := M[1, 1] * M[3, 1] + M[1, 2] * M[3, 2] + M[1, 3] * M[3, 3] + M[1, 4] * M[3, 4];
S14 := M[1, 1] * M[4, 1] + M[1, 2] * M[4, 2] + M[1, 3] * M[4, 3] + M[1, 4] * M[4, 4];
S23 := M[2, 1] * M[3, 1] - M[2, 2] * M[3, 2] - M[2, 3] * M[3, 3] - M[2, 4] * M[3, 4];
S24 := M[2, 1] * M[4, 1] - M[2, 2] * M[4, 2] - M[2, 3] * M[4, 3] - M[2, 4] * M[4, 4];
S34 := M[3, 1] * M[4, 1] - M[3, 2] * M[4, 2] - M[3, 3] * M[4, 3] - M[3, 4] * M[4, 4];
S := expand(-S12^2 - S13^2 - S14^2 + S23^2 + S24^2 + S34^2);
tr4_9 := expand(tr4_8 + 16 * S);

G_tail_total := add(G_tail);
L := expand(G_tail_total * squares1x);
R:=expand((M[1, 1]^2 - squaresx1)^2);
tr4_10 := tr4_9 + 16 * L - 8 * R;
tr4_11 := tr4_10 + 16 * M[1, 1]^4 + 8 * (M[2, 1]^4 + M[3, 1]^4 + M[4, 1]^4);
