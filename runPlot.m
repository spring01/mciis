
cart = [...
    6                   0.000000    0.463133    0.000000
    6                  -0.935411   -0.720058    0.000000
    8                   1.208199    0.386722    0.000000
    1                  -0.503791    1.457138    0.000000
    1                  -0.369919   -1.654663    0.000000
    1                  -1.589710   -0.677352    0.880800
    1                  -1.589710   -0.677352   -0.880800];

mol = Molecule(cart);
basisSet = '6-31g*';
dft = 'b3lypv5';

matpsi = MatPsi2(mol.cartesian, basisSet, 0, 1);
scf = RKS(RHF.MatPsi2Interface(matpsi), dft);


load aldeOut.mat;

% load sih3hOut.mat
% output.energySet{1} = output.energySet{1}(1:40);

% load cdimOut.mat;
% output.energySet{1} = output.energySet{1}(1:40);
% output.energySet{2} = output.energySet{2}(1:40);


scf.Plot(output);




