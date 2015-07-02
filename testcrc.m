
cart = [...
 24                0.00000000    0.00000000    0.40000000
 6                 0.00000000    0.00000000   -1.60000000];

mol = Molecule(cart);
basisSet = '6-31g';
dft = 'b3lyp';

matpsi = MatPsi2(mol.cartesian, basisSet, 0, 3);
% scf = UHF(RHF.MatPsi2Interface(matpsi));
scf = UKS(RHF.MatPsi2Interface(matpsi), dft);
output = scf.Run();
scf.Plot(output);

