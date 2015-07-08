cart = [...
    92                  0.000000    0.000000    0.000000
    9                   1.143154    1.143154    1.143154
    9                  -1.143154   -1.143154    1.143154
    9                   1.143154   -1.143154   -1.143154
    9                  -1.143154    1.143154   -1.143154];

info.chargeMult = [0 3];
info.cartesian = cart;
info.method = 'b3lyp';
info.basisSet = 'genecp';
info.ecpFile = 'uf4.ecp';
scf = G09USCF(info);
output = scf.Run();
scf.Plot(output);


% allTypes = {'C20', 'M20', 'ECe20', 'EMe20', 'ACe20', 'AMe20'};

