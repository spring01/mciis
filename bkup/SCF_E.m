function [energy, iter] = SCF_E(obj, guessDensity)
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);
densVec = guessDensity;

% diis
ediis = obj.EDIIS(20);
energy = 0;
fockVec = [];
for iter = 1:obj.maxSCFIter
    oldEnergy = energy;
    oldFockVec = fockVec;
    fockVec = obj.DensVecToFockVec(densVec);
    energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    
    % damping at the 2nd iteration
    if(iter == 2)
        dampingCoeff = 0.2;
        fockVec = obj.Damping(dampingCoeff, fockVec, oldFockVec);
        densVec = obj.Damping(dampingCoeff, densVec, oldDensVec);
        energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    end
    
    % diis extrapolate Fock matrix
    ediis.Push(fockVec, densVec, energy);
    fockVec = ediis.OptFockVector();
    disp('ediis')
    
%     disp(energy);
    
    oldDensVec = densVec;
    [densVec, orbital] = obj.SolveFockVec(fockVec, inv_S_Half);
    
    if(mean(sqrt(mean((densVec - oldDensVec).^2))) < obj.RMSDensityThreshold ...
            && mean(max(abs(densVec - oldDensVec))) < obj.MaxDensityThreshold ...
            && abs(energy - oldEnergy) < obj.EnergyThreshold)
        break;
    end
    
    disp(['done iter ', num2str(iter)])
end

obj.orbital = orbital;
obj.densVec = densVec;

end




