info.chargeMult = [0 1];
info.cartesian = [...
    8 0.0 0.0 0.0
    1 0.0 0.0 1.0
    1 0.0 1.0 0.0];
info.method = 'hf';
info.basisSet = '6-31g*';

scf = G09RSCF(info);
[guessDensity, guessOrbital] = scf.CoreGuess();

[ener1, energySet1, iter1] = scf.SCF(guessOrbital, 'C20');

