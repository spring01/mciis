cart = [...
    92                  0.000000    0.000000    0.000000
    9                   1.143154    1.143154    1.143154
    9                  -1.143154   -1.143154    1.143154
    9                   1.143154   -1.143154   -1.143154
    9                  -1.143154    1.143154   -1.143154];

% cart = [...
%     14                  0.000000    0.000000    0.140556
%     1                   0.000000    1.385929    0.630556
%     1                   1.200250   -0.692965    0.630556
%     1                  -1.200250   -0.692965    0.630556
%     1                   0.000000    0.000000   -3.859444];

info.chargeMult = [0 1];
info.cartesian = cart;
info.method = 'b3lyp';
info.basisSet = 'genecp';
info.ecpFile = 'uf4.ecp';

scf = G09RSCF(info);
[guessDensity, guessOrbital] = scf.CoreGuess();

[ener1, energySet1, iter1] = scf.SCF(guessOrbital, 'LISTd5');

fprintf('%0.8f  %d \n',ener1, iter1);