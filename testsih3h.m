
cart = [...
    14                  0.000000    0.000000    0.140556
    1                   0.000000    1.385929    0.630556
    1                   1.200250   -0.692965    0.630556
    1                  -1.200250   -0.692965    0.630556
    1                   0.000000    0.000000   -3.859444];

% cart = [...
%     14                -1.55096010   -0.68685376    0.00000000
%     1                 -1.06097692   -2.07278899    0.00000000
%     1                 -1.06095162    0.00610450    1.20025020
%     1                 -1.06095162    0.00610450   -1.20025020
%     1                 -5.55096010   -0.68680447    0.00000000];



mol = Molecule(cart);
basisSet = '6-31g*';
dft = 'svwn5';

matpsi = MatPsi2(mol.cartesian, basisSet, 0, 1);
scf = RHF(RHF.MatPsi2Interface(matpsi));
% scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
output = scf.Run();
scf.Plot(output);

