classdef UKS < UHF & RKS
    
    methods
        
        function obj = UKS(properties, dft)
            obj@RKS(properties, dft);
            obj@UHF(properties);
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            occOrbAlpha = orbital{1}(:, 1:obj.numElectrons(1));
            occOrbBeta = orbital{2}(:, 1:obj.numElectrons(2));
            jMat = obj.matpsi2.JK_OccOrbToJ(occOrbAlpha, occOrbBeta);
            obj.currentV = obj.matpsi2.DFT_OccOrbToV(occOrbAlpha, occOrbBeta);
            gMat = repmat(sum(jMat, 3), [1 1 2]) + obj.currentV;
            if(obj.hfExcCoeff ~= 0)
                kMat = obj.matpsi2.JK_OccOrbToK(occOrbAlpha, occOrbBeta);
                gMat = gMat - obj.hfExcCoeff * kMat;
            end
            oeiVec = reshape(obj.coreHamilt, [], 1);
            fockVec = repmat(oeiVec, 1, 2) + ...
                [reshape(gMat(:,:,1), [], 1), reshape(gMat(:,:,2), [], 1)];
        end
        
        function elecEnergy = ElecEnergy(obj, fockVec, densVec)
            elecEnergy = obj.ElecEnergy@UHF(fockVec, densVec) ...
                - reshape(obj.currentV(:,:,1), 1, []) * densVec(:, 1) / 2 ...
                - reshape(obj.currentV(:,:,2), 1, []) * densVec(:, 2) / 2 ...
                + obj.matpsi2.DFT_EnergyXC();
        end
        
    end
    
end