% cart = [...
%  8                 -2.54062035   -0.22895125    0.00000000
%  1                 -1.58062035   -0.22895125    0.00000000
%  1                 -2.86107494    0.67598458    0.00000000];

cart = [...
    8 0.0 0.0 0.0
    1 0.0 0.0 1.0
    1 0.0 1.0 0.0];

% cart = [...
%  14                -1.55096010   -0.68685376    0.00000000
%  1                 -1.06097692   -2.07278899    0.00000000
%  1                 -1.06095162    0.00610450    1.20025020
%  1                 -1.06095162    0.00610450   -1.20025020
%  1                 -5.55096010   -0.68680447    0.00000000];

% cart = [...
%     14                  0.000000    0.000000    0.140556
%     1                   0.000000    1.385929    0.630556
%     1                   1.200250   -0.692965    0.630556
%     1                  -1.200250   -0.692965    0.630556
%     1                   0.000000    0.000000   -3.859444];

% cart = [...
%  24                0.00000000    0.00000000    0.40000000
%  6                 0.00000000    0.00000000   -1.60000000];


% cart = [...
%     6                 2.789343    0.740153   -0.000006
%     1                 3.699705    1.315321   -0.000005
%     6                 1.489141    1.148729    0.000011
%     1                 1.095680    2.150351    0.000035
%     6                 1.505673   -1.105101   -0.000008
%     1                 1.220972   -2.143783    0.000073
%     7                 2.771491   -0.656327   -0.000009
%     1                 3.603538   -1.250416    0.000049
%     48               -1.426682   -0.000879    0.000000
%     7                 0.679350   -0.019671   -0.000010];

% cart = [...
% 6                  -1.171138   -0.148546    0.000004
% 1                  -1.714145    0.219764   -0.880470
% 1                  -1.713691    0.218958    0.881075
% 1                  -1.155613   -1.240682   -0.000601
% 6                   0.234078    0.399597    0.000022
% 1                   0.301980    1.511628   -0.000146
% 8                   1.237979   -0.276996   -0.000002];

mol = Molecule(cart);
basisSet = '6-31g*';
dft = 'b3lyp';
diisType = 'C20';

matpsi = MatPsi2(mol.cartesian, basisSet, 0, 3);
% scf = UHF(RHF.MatPsi2Interface(matpsi));
scf = UKS(RHF.MatPsi2Interface(matpsi), dft);
[guessDensity, guessOrbital] = scf.CoreGuess();

[ener1, energySet1, iter1] = scf.SCF(guessOrbital, diisType);

fprintf('%0.8f  %d \n',ener1, iter1);

figure();
hold();
plot(log10(abs(energySet1 - ener1)), 'r');
scatter(1:length(energySet1), log10(abs(energySet1 - ener1)), 72, 'square', 'r', 'filled');




