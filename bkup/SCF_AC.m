function [energy, iter] = SCF_AC(obj, guessDensity)
nbf = size(obj.overlapMat, 1);

oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = reshape(guessDensity, [], 1);

% diis
adiis = ADIIS(obj.coreHamilt, 20);
cdiis = CDIIS(obj.overlapMat, 20);
energy = 0;
maxErrSet = [];
fockVec = zeros(size(oeiVec));
for iter = 1:obj.maxSCFIter
    oldEnergy = energy;
    oldFockVec = fockVec;
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    
    % damping at the 2nd iteration
    if(iter == 2)
        coeffs = [0.25; 0.75];
        fockVec = [fockVec, oldFockVec] * coeffs;
        densVec = [densVec, oldDensVec] * coeffs;
        energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    end
    
    % diis extrapolate Fock matrix
    adiis.Push(fockVec, densVec);
    cdiis.Push(fockVec, densVec);
    maxErr = cdiis.MaxError();
%     disp(maxErr);
    if(maxErr ~= 0)
        maxErrSet(iter) = maxErr;
    else
        maxErrSet(iter) = 1;
    end
    
    if(maxErr > 1e-4)
        [fockVecSet, coeffs] = adiis.FockCoeffs();
        disp('adiis')
    else
        [fockVecSet, coeffs] = cdiis.FockCoeffs();
        disp('cdiis')
    end
    fockVec = fockVecSet * coeffs;
    
%     disp(energy);
    
    oldDensVec = densVec;
    [densVec, orbital] ...
        = obj.DiagonalizeFock(reshape(fockVec, nbf, []), ...
        inv_S_Half);
    
    if(sqrt(mean((densVec - oldDensVec).^2)) < obj.RMSDensityThreshold ...
            && max(abs(densVec - oldDensVec)) < obj.MaxDensityThreshold ...
            && abs(energy - oldEnergy) < obj.EnergyThreshold)
        break;
    end
    
    disp(['done iter ', num2str(iter)])
end

obj.orbital = orbital;
obj.densVec = densVec;

end




