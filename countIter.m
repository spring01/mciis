load gdb11randpertstd.mat

basisSet = '6-31+g*';
dft = 'b3lyp';
diisType = 'EM20';

totalNum = length(molecules);

ener = zeros(totalNum, 1);
iter = zeros(totalNum, 1);
for ind = 25:25
    cart = molecules{ind}.cartesian;
    
    mol = Molecule(cart);
    
    matpsi = MatPsi2(mol.cartesian, basisSet, 0, 1);
    matpsi.SCF_SetSCFType('rks');
    matpsi.JK_Initialize('dfjk');
    scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
    
    
    [guessDensity, guessOrbital] = scf.CoreGuess();
    [ener(ind), energySet1, iter(ind)] = scf.SCF(guessOrbital, diisType);
    disp(ind);
end

output.ener = ener;
output.iter = iter;
% save('pertM20.mat', 'output');



