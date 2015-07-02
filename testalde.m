

cart = [...
6       0.0000      0.4631      0.0000
6      -0.9354     -0.7201      0.0000
8       1.2082      0.3867      0.0000
1      -0.5038      1.4571      0.0000
1      -0.3699     -1.6547      0.0000
1      -1.5897     -0.6774      0.8808
1      -1.5897     -0.6774     -0.8808];

mol = Molecule(cart);
basisSet = '6-31g*';
dft = 'b3lypv5';

matpsi = MatPsi2(mol.cartesian, basisSet, 0, 1);
scf = RHF(RHF.MatPsi2Interface(matpsi));
% scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
output = scf.Run();
scf.Plot(output);

