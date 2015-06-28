function [energy, iter] = SCF_ECMmix(obj, guessDensity)
nbf = size(obj.overlapMat, 1);

oeiVec = reshape(obj.coreHamilt, [], 1);
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);

densVec = reshape(guessDensity, [], 1);

% diis
ediis = EDIIS(obj.coreHamilt, 20);
cdiis = CDIIS(obj.overlapMat, 20);
mciis = MCIIS(obj.overlapMat, 20);
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
    ediis.Push(fockVec, densVec, energy);
    cdiis.Push(fockVec, densVec);
    mciis.Push(fockVec, densVec);
    maxErr = cdiis.MaxError();
%     disp(maxErr);
    if(maxErr ~= 0)
        maxErrSet(iter) = maxErr;
    else
        maxErrSet(iter) = 1;
    end
    if(maxErr > 1e-1 || maxErr > 1.1 * min(maxErrSet))
        [fockVecSet, coeffs] = ediis.FockCoeffs();
        disp('ediis');
    elseif(maxErr < 1e-4)
        [fockVecSet, coeffs] = mciis.FockCoeffs();
        disp('mciis');
    else
        [~, ediisCoeffs] = ediis.FockCoeffs();
        [fockVecSet, cdiisCoeffs] = cdiis.FockCoeffs();
        coeffs = 10.*maxErr .* ediisCoeffs + (1 - 10.*maxErr) .* cdiisCoeffs;
        disp('mix cdiis ediis');
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




