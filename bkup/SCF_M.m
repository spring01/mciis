function [energy, iter] = SCF_M(obj, guessDensity)
nbf = size(obj.overlapMat, 1);

oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = reshape(guessDensity, [], 1);

% diis
mciis = MCIIS(obj.overlapMat, 20);
energy = 0;
fockVec = oeiVec;
for iter = 1:obj.maxSCFIter
    oldFockVec = fockVec;
    oldEnergy = energy;
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
    mciis.Push(fockVec, densVec);
    fockVec = mciis.Extrapolate();
    disp('mciis');
    
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




