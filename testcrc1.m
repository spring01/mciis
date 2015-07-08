
cart = [...
 24                0.00000000    0.00000000    0.40000000
 6                 0.00000000    0.00000000   -1.60000000];

mol = Molecule(cart);
basisSet = '6-31g';
dft = 'b3lyp';

matpsi = MatPsi2(mol.cartesian, basisSet, 0, 1);
scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
output = scf.Run();
scf.Plot(output);


% allTypes = {'C20', 'M20', 'ECe20', 'EMe20', 'ACe20', 'AMe20'};


