

cart = [...
    48     0.000000     0.000000     0.000000
    7      0.000000     0.000000    -2.260001
    7     -0.685444     0.000000    -4.348035
    6      0.676053     0.000000    -4.385069
    6      1.085240     0.000000    -3.091231
    6     -1.044752     0.000000    -3.060220
    1      1.231530     0.000000    -5.300759
    1      2.088641     0.000000    -2.711077
    1     -2.068750     0.000000    -2.726515
    1     -1.313170     0.000000    -5.174718];

mol = Molecule(cart);
basisSet = '3-21g';
dft = 'b3lyp';

matpsi = MatPsi2(mol.cartesian, basisSet, 2, 1);
% scf = RHF(RHF.MatPsi2Interface(matpsi));
scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
output = scf.Run();
scf.Plot(output);

