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
corner_det := Determinant([[M[2, 2], M[2, 3], M[2, 4]], [M[3, 2], M[3, 3], M[3, 4]], [M[4, 2], M[4, 3], M[4, 4]]]); # 9 mult

print("Verification of the formula for tr(H^3): ", expand(tr3 - 3 * M[1, 1] * tr2 + 8 * M[1, 1]^3 - 24 * (M[2, 1] * F2 + M[3, 1] * F3 + M[4, 1] * F4 + corner_det))); # 5 mult

# for tr3 --> 23 mult

full_det := Determinant(Matrix([seq([seq(M[i, j], j=1..4)], i=1..4)]));

tr4_1 := expand(tr4 + 32 * full_det);

tr4_2 := expand(tr4_1 - 12 * (tr2 / 4)^2);

G := [seq(add([seq(M[i, j]^2, j=1..4)]), i=1..4)];
G2 := expand(add([seq(G[i]^2, i=1..4)]));

G_tail := [seq(add([seq(M[i, j]^2, j=2..4)]), i=2..4)];
G_tail2 := expand(add([seq(G_tail[i]^2, i=1..3)]));

tr4_3 := expand(tr4_2 + 8 * G_tail2);

G_tail_op := [seq(add([seq(M[j, i]^2, j=2..4)]), i=2..4)];
G_tail_op2 := expand(add([seq(G_tail_op[i]^2, i=1..3)]));
sum4 := add(seq(M[i, j]^4, i=2..4, j=2..4));

tr4_4 := tr4_3 + 8 * G_tail_op2 - 8 * sum4;

anticommuting_with1 := expand(M[1, 2]^2 * (G_tail_op[2] + G_tail_op[3]) + M[1, 3]^2 * (G_tail_op[1] + G_tail_op[3]) + M[1, 4]^2 * (G_tail_op[1] + G_tail_op[2]));
anticommuting_with1_op := expand(M[2, 1]^2 * (G_tail[2] + G_tail[3]) + M[3, 1]^2 * (G_tail[1] + G_tail[3]) + M[4, 1]^2 * (G_tail[1] + G_tail[2]));

tr4_5 := tr4_4 + 16 * anticommuting_with1 + 16 * anticommuting_with1_op;

squares1x := add([seq(M[1, i]^2, i=2..4)]);
squaresx1 := add([seq(M[i, 1]^2, i=2..4)]);

tr4_6 := tr4_5 + 8 * expand(squares1x^2 + squaresx1^2);

tr4_7 := tr4_6 - 128 * expand(M[1, 1] * corner_det);
