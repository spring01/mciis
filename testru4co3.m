cart = [...
    44                  -0.219401   -1.353713   -0.552094
    44                  -0.218751    1.354097   -0.551737
    44                   1.218036   -0.000363    0.874787
    44                  -1.831167    0.000061    0.392072
    6                    2.757624   -0.000445   -0.123186
    8                    3.713844   -0.000111   -0.804262];

info.chargeMult = [0 3];
info.cartesian = cart;
info.method = 'b3lyp';
info.basisSet = 'lanl2dz';
scf = G09USCF(info);
output = scf.Run();
scf.Plot(output);


% allTypes = {'ECe20', 'EMe20', 'ACe20', 'AMe20'};

