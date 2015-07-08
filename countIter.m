load gdb13rand100g09std.mat

basisSet = '6-31+g*';
dft = 'b3lyp';
diisType = 'C20';

totalNum = length(molecules);

ener = zeros(totalNum, 1);
iter = zeros(totalNum, 1);
for ind = 1:totalNum
    cart = molecules{ind}.cartesian;
    
    mol = Molecule(cart);
    
    if(mod(sum(cart(:, 1)), 2) == 0)
        matpsi = MatPsi2(mol.cartesian, basisSet, 0, 1);
        matpsi.SCF_SetSCFType('rks');
        scf = RKS(RHF.MatPsi2Interface(matpsi), dft);
    else
        matpsi = MatPsi2(mol.cartesian, basisSet, 0, 2);
        matpsi.SCF_SetSCFType('uks');
        scf = UKS(RHF.MatPsi2Interface(matpsi), dft);
    end
    
    
    [guessDensity, guessOrbital] = scf.CoreGuess();
    [ener(ind), energySet1, iter(ind)] = scf.SCF(guessOrbital, diisType);
    disp(ind);
end

output.ener = ener;
output.iter = iter;
save('nopert.mat', 'output');



