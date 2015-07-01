classdef RKS < RHF
    
    properties (Access = protected)
        
        currentV;
        hfExcCoeff;
        
    end
    
    methods
        
        function obj = RKS(properties, dft)
            obj = obj@RHF(properties);
            obj.matpsi2.DFT_Initialize(dft);
            if(strcmpi(dft, 'b3lyp') || strcmpi(dft, 'b3lypv5'))
                obj.hfExcCoeff = 0.2;
            else
                obj.hfExcCoeff = 0;
            end
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            occOrb = orbital(:, 1:obj.numElectrons/2);
            obj.currentV = obj.matpsi2.DFT_OccOrbToV(occOrb);
            gMat = 2 .* obj.matpsi2.JK_OccOrbToJ(occOrb) + obj.currentV;
            if(obj.hfExcCoeff ~= 0)
                gMat = gMat - obj.hfExcCoeff * obj.matpsi2.JK_OccOrbToK(occOrb);
            end
            fockVec = reshape(obj.coreHamilt, [], 1) + reshape(gMat, [], 1);
        end
        
        function elecEnergy = SCFEnergy(obj, fockVec, densVec)
            elecEnergy = obj.SCFEnergy@RHF(fockVec, densVec) ...
                - reshape(obj.currentV, 1, []) * densVec ...
                + obj.matpsi2.DFT_EnergyXC();
        end
        
    end
    
end