function [energy, iter] = SCF_ECmix(obj, guessDensity)

oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = reshape(guessDensity, [], 1);

% diis
ediis = obj.EDIIS(20);
cdiis = obj.CDIIS(20);
energy = 0;
maxErrSet = [];
fockVec = zeros(size(oeiVec));
for iter = 1:obj.maxSCFIter
    oldEnergy = energy;
    oldFockVec = fockVec;
    fockVec = oeiVec + obj.DensVecToGVec(densVec);
    energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    
    % damping at the 2nd iteration
    if(iter == 2)
        dampingCoeff = 0.20;
        fockVec = obj.Damping(dampingCoeff, fockVec, oldFockVec);
        densVec = obj.Damping(dampingCoeff, densVec, oldDensVec);
        energy = obj.ElecEnergy(fockVec, densVec) + obj.nucRepEnergy;
    end
    
    % diis extrapolate Fock matrix
    obj.EDIISPush(ediis, fockVec, densVec, energy);
    obj.CDIISPush(cdiis, fockVec, densVec);
    maxErr = obj.CDIISMaxError(cdiis);
%     disp(maxErr);
    if(maxErr ~= 0)
        maxErrSet(iter) = maxErr;
    else
        maxErrSet(iter) = 1;
    end
    
    if(maxErr > 1e-1 || maxErr > 1.1 * min(maxErrSet))
        [fockVecSet, coeffs] = obj.DIISFockCoeffs(ediis);
        disp('ediis')
    elseif(maxErr < 1e-4)
        [fockVecSet, coeffs] = obj.DIISFockCoeffs(cdiis);
        disp('cdiis')
    else
        [~, ediisCoeffs] = obj.DIISFockCoeffs(ediis);
        [fockVecSet, cdiisCoeffs] = obj.DIISFockCoeffs(cdiis);
        coeffs = 10.*maxErr .* ediisCoeffs + (1 - 10.*maxErr) .* cdiisCoeffs;
        disp('mix cdiis ediis');
    end
    fockVec = obj.DIISCalcFockVec(fockVecSet, coeffs);
    
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




