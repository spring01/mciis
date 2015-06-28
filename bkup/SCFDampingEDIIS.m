function [energy, iter] = SCFDampingEDIIS(obj)
nbf = size(obj.overlapMat, 1);

oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = obj.DiagonalizeFock(reshape(oeiVec, nbf, []), inv_S_Half);

% diis
ediis = EDIIS(obj.coreHamilt, 20);
energy = 0;
fockVec = oeiVec;
for iter = 1:obj.maxSCFIter
    oldFockVec = fockVec;
    oldEnergy = energy;
    fockVec = oeiVec + reshape(obj.DensToG(reshape(densVec, nbf, [])), [], 1);
    energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    
    % diis extrapolate Fock matrix
    ediis.Push(fockVec, densVec, energy);
    if(iter == 2)
        fockVec = 0.25 * fockVec + 0.75 * oldFockVec;
        densVec = 0.25 * densVec + 0.75 * oldDensVec;
        energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    else
        fockVec = ediis.Interpolate();
    end
    
    disp(energy);
    
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




