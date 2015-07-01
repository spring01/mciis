info.chargeMult = [0 1];
info.cartesian = [...
    8 0.0 0.0 0.0
    1 0.0 0.0 1.0
    1 0.0 1.0 0.0];
info.method = 'b3lyp';
info.basisSet = 'sto-3g';

G09RSCF(info)