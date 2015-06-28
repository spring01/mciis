function [energy, iter] = SCF_ECe(obj, guessDensity)
nbf = size(obj.overlapMat, 1);

oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = reshape(guessDensity, [], 1);

% diis
ediis = EDIIS(obj.coreHamilt, 20);
cdiis = CDIIS(obj.overlapMat, 20);
energy = 0;
energyDiffSet = [];
fockVec = zeros(size(oeiVec));
for iter = 1:obj.maxSCFIter
    oldEnergy = energy;
    oldFockVec = fockVec;
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    
    % damping at the 2nd iteration
    if(iter == 2)
        dampingCoeff = 0.25;
        coeffs = [dampingCoeff; (1 - dampingCoeff)];
        fockVec = [fockVec, oldFockVec] * coeffs;
        densVec = [densVec, oldDensVec] * coeffs;
        energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    end
    
    % diis extrapolate Fock matrix
    ediis.Push(fockVec, densVec, energy);
    cdiis.Push(fockVec, densVec);
    energyDiff = abs(energy - oldEnergy);
%     disp(maxErr);
    if(energyDiff ~= 0)
        energyDiffSet(iter) = energyDiff;
    else
        energyDiffSet(iter) = 1;
    end
    
%     if(energyDiff > 1e-3 || energyDiff > 1.1 * min(energyDiffSet))
    if(energyDiff > 1e-1)
        [fockVecSet, coeffs] = ediis.FockCoeffs();
        disp('ediis')
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




