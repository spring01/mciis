
cart = [...
    14                  0.000000    0.000000    0.140556
    1                   0.000000    1.385929    0.630556
    1                   1.200250   -0.692965    0.630556
    1                  -1.200250   -0.692965    0.630556
    1                   0.000000    0.000000   -3.859444];



mol = Molecule(cart);
basisSet = '6-31g*';
dft = 'svwn5';

matpsi = MatPsi2(mol.cartesian, basisSet, 0, 1);
scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
output = scf.Run();
scf.Plot(output);


% allTypes = {'C6', 'C20', 'M6', 'M20', 'ECe6', 'ECe20', 'EMe6', 'EMe20'};

