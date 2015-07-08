
load gdb7g09stdgeom.mat

% randNum = ceil(rand() * length(molecules))
randNum = 2
% randNum = 90;

% randNum = 910;
cart = molecules{randNum}.cartesian;
randMat = 2*(rand(size(cart(:, 2:end))) - 0.5);
cart(:, 2:end) = cart(:, 2:end) + randMat;

mol = Molecule(cart);
basisSet = '6-31g*';
dft = 'b3lyp';

if(mod(mol.NumElectrons(), 2))
    matpsi = MatPsi2(mol.cartesian, basisSet, 0, 2);
    matpsi.SCF_SetSCFType('uks');
    scf = UKS(RHF.MatPsi2Interface(matpsi), dft);
else
    matpsi = MatPsi2(mol.cartesian, basisSet, 0, 1);
    matpsi.SCF_SetSCFType('rks');
    scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
end

% matpsi = MatPsi2(mol.cartesian, basisSet, 0, mult);
% matpsi.SCF_SetSCFType('uks');
% scf = RHF(RHF.MatPsi2Interface(matpsi));
% scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
output = scf.Run();
scf.Plot(output);

