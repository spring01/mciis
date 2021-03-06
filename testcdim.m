

cart = [...
    48                 -0.019657   -1.492605    0.000000
    7                   0.000000    0.767310    0.000000
    7                   0.703579    2.849304    0.000000
    6                  -0.657544    2.898178    0.000000
    6                  -1.077969    1.607948    0.000000
    6                   1.051673    1.558412    0.000000
    1                  -1.205035    3.818665    0.000000
    1                  -2.084639    1.236536    0.000000
    1                   2.072729    1.215813    0.000000
    1                   1.338472    3.670495    0.000000];

mol = Molecule(cart);
basisSet = '3-21g';
dft = 'b3lyp';

matpsi = MatPsi2(mol.cartesian, basisSet, 2, 1);
scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
output = scf.Run();
scf.Plot(output);


% allTypes = {'C6', 'C20', 'M6', 'M20', 'ECe6', 'ECe20', 'EMe6', 'EMe20'};


