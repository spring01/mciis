function [energy, energySet, iter] = SCF(obj, guessOrbital, diisType)
inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);
orbital = guessOrbital;
densVec = obj.OrbToDensVec(orbital);

% diis
cdiis20 = obj.CDIIS(20);
ediis20 = obj.EDIIS(20);
mciis20 = obj.MCIIS(20);

% initialize some variables
energy = 0;
fockVec = [];
maxErrSet = [];
energySet = [];
for iter = 1:obj.maxSCFIter
    oldFockVec = fockVec;
    oldEnergy = energy;
    fockVec = obj.OrbToFockVec(orbital);
    energy = obj.SCFEnergy(fockVec, densVec);
    
%     % damping at the 2nd iteration
%     if(iter == 2)
%         dampingCoeff = 0.25;
%         fockVec = obj.Damping(dampingCoeff, fockVec, oldFockVec);
%         densVec = obj.Damping(dampingCoeff, densVec, oldDensVec);
%         energy = obj.SCFEnergy(fockVec, densVec);
%     end
    
    % diis extrapolate Fock matrix
    cdiis20.Push(fockVec, densVec);
    ediis20.Push(fockVec, densVec, energy);
    mciis20.Push(fockVec, densVec);
    switch(diisType)
        case('ECmix20')
            [fockVec, maxErrSet] = ECmix(ediis20, cdiis20, maxErrSet, iter);
        case('C20')
            fockVec = cdiis20.OptFockVector();
            disp('cdiis(20)');
        case('E20')
            fockVec = ediis20.OptFockVector();
            disp('ediis(20)')
        case('EC20')
            [fockVec, maxErrSet] = EC(ediis20, cdiis20, maxErrSet, iter);
        case('ECe20')
            fockVec = ECe(ediis20, cdiis20, abs(energy - oldEnergy));
        case('M20')
            fockVec = mciis20.OptFockVector();
            disp('mciis(20)');
        case('EM20')
            [fockVec, maxErrSet] = EM(ediis20, cdiis20, mciis20, maxErrSet, iter);
        case('EMe20')
            fockVec = EMe(ediis20, mciis20, abs(energy - oldEnergy));
    end
    
    energySet(iter) = energy; %#ok
%     disp(energy);
    
    oldDensVec = densVec;
    orbital = obj.SolveFockVec(fockVec, inv_S_Half);
    densVec = obj.OrbToDensVec(orbital);
    
    if(mean(sqrt(mean((densVec - oldDensVec).^2))) < obj.RMSDensityThreshold ...
            && mean(max(abs(densVec - oldDensVec))) < obj.MaxDensityThreshold ...
            && abs(energy - oldEnergy) < obj.EnergyThreshold)
        break;
    end
    
    disp(['done iter ', num2str(iter)])
end

end


function [fockVec, maxErrSet] = ECmix(ediis, cdiis, maxErrSet, iter)
maxErr = cdiis.MaxError();
if(maxErr ~= 0)
    maxErrSet(iter) = maxErr;
else
    maxErrSet(iter) = 1;
end
if(maxErr > 1e-1 || maxErr > 1.1 * min(maxErrSet))
    fockVec = ediis.OptFockVector();
    disp('ediis')
elseif(maxErr < 1e-4)
    fockVec = cdiis.OptFockVector();
    disp('cdiis')
else
    [~, ediisCoeffs, ~] = ediis.OptFockVector();
    [~, cdiisCoeffs, fockVecSet] = cdiis.OptFockVector();
    coeffs = 10.*maxErr .* ediisCoeffs + (1 - 10.*maxErr) .* cdiisCoeffs;
    fockVec = cdiis.CalcFockVec(coeffs, fockVecSet);
    disp('mix cdiis ediis');
end
end

function [fockVec, maxErrSet] = EC(ediis, cdiis, maxErrSet, iter)
maxErr = cdiis.MaxError();
if(maxErr ~= 0)
    maxErrSet(iter) = maxErr;
else
    maxErrSet(iter) = 1;
end
% if(maxErr > 1e-1 || maxErr > 1.1 * min(maxErrSet))
if(maxErr > 1e-2)
    fockVec = ediis.OptFockVector();
    disp('ediis')
else
    fockVec = cdiis.OptFockVector();
    disp('cdiis')
end
end

function fockVec = ECe(ediis, cdiis, energyDiff)
if(energyDiff > 1e-2)
    fockVec = ediis.OptFockVector();
    disp('ediis')
else
    fockVec = cdiis.OptFockVector();
    disp('cdiis')
end
end

function [fockVec, maxErrSet] = EM(ediis, cdiis, mciis, maxErrSet, iter)
maxErr = cdiis.MaxError();
if(maxErr ~= 0)
    maxErrSet(iter) = maxErr;
else
    maxErrSet(iter) = 1;
end
% if(maxErr > 1e-1 || maxErr > 1.1 * min(maxErrSet))
if(maxErr > 1e-2)
    fockVec = ediis.OptFockVector();
    disp('ediis')
else
    fockVec = mciis.OptFockVector();
    disp('mciis')
end
end

function fockVec = EMe(ediis, mciis, energyDiff)
if(energyDiff > 1e-2)
    fockVec = ediis.OptFockVector();
    disp('ediis')
else
    fockVec = mciis.OptFockVector();
    disp('mciis')
end
end


