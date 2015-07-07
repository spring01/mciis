% cart = [...
%     92                  0.000000    0.000000    0.000000
%     9                   1.143154    1.143154    1.143154
%     9                  -1.143154   -1.143154    1.143154
%     9                   1.143154   -1.143154   -1.143154
%     9                  -1.143154    1.143154   -1.143154];

% cart = [...
%     14                  0.000000    0.000000    0.140556
%     1                   0.000000    1.385929    0.630556
%     1                   1.200250   -0.692965    0.630556
%     1                  -1.200250   -0.692965    0.630556
%     1                   0.000000    0.000000   -3.859444];

cart = [...
    44                  -0.219401   -1.353713   -0.552094
    44                  -0.218751    1.354097   -0.551737
    44                   1.218036   -0.000363    0.874787
    44                  -1.831167    0.000061    0.392072
    6                    2.757624   -0.000445   -0.123186
    8                    3.713844   -0.000111   -0.804262];

info.chargeMult = [0 1];
info.cartesian = cart;
info.method = 'b3lyp';
info.basisSet = 'lanl2dz';
% info.ecpFile = 'uf4.ecp';
scf = G09RSCF(info);
[guessDensity, guessOrbital] = scf.CoreGuess();
[ener1, energySet1, iter1, orbital] = scf.SCF(guessOrbital, 'ECe20');

% EC20 -488.70476735  60 
% EM20 -488.70476735  63 

% ECe20 -488.70476735  67 
% EMe20 -488.70476735  61 

fprintf('%0.8f  %d \n',ener1, iter1);

figure();
hold();
plot(log10(abs(energySet1 - ener1)), 'r');
scatter(1:length(energySet1), log10(abs(energySet1 - ener1)), 72, 'square', 'r', 'filled');

