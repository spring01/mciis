chargeMult = [0 1];
cart = [...
    8 0.0 0.0 0.0
    1 0.0 0.0 1.0
    1 0.0 1.0 0.0];
method = 'b3lyp';
basisSet = 'sto-3g';
orbital = rand(7) - 0.5;

fileStr = G09InputStr(chargeMult, cart, method, basisSet, orbital);
gjfFile = fopen('temp.gjf', 'w');
fprintf(gjfFile, '%s', fileStr);
fclose(gjfFile);
system('g09 temp.gjf');

G09FileIsValid()

% G09ReadMatrix({'overlap', 'coreHamilt', 'fockAlpha'})
