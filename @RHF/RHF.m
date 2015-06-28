classdef RHF < handle
    
    properties (SetAccess = protected)
        
        energySet;
        
    end
    
    properties (Access = protected)
        
        overlapMat;
        coreHamilt;
        nucRepEnergy;
        numElectrons;
        matpsi2;
        
        maxSCFIter = 100;
        RMSDensityThreshold = 1e-8;
        MaxDensityThreshold = 1e-6;
        EnergyThreshold = 1e-6;
        
    end
    
    methods
        
        function obj = RHF(properties)
            obj.overlapMat = properties.overlapMat;
            obj.coreHamilt = properties.coreHamilt;
            obj.nucRepEnergy = properties.nucRepEnergy;
            obj.numElectrons = properties.numElectrons;
            obj.matpsi2 = properties.matpsi2;
        end
        
        function [densVec, orbital] = CoreGuess(obj)
            coreFockVec = reshape(obj.coreHamilt, [], 1);
            inv_S_Half = eye(size(obj.overlapMat)) / sqrtm(obj.overlapMat);
            orbital = obj.SolveFockVec(coreFockVec, inv_S_Half);
            densVec = obj.OrbToDensVec(orbital);
        end
        
    end
    
    methods (Access = protected)
        
        function densVec = OrbToDensVec(obj, orbital)
            occOrb = orbital(:, 1:obj.numElectrons/2);
            densVec = reshape(occOrb * occOrb', [], 1);
        end
        
        function fockVec = OrbToFockVec(obj, orbital)
            occOrb = orbital(:, 1:obj.numElectrons/2);
            gMat = 2 .* obj.matpsi2.JK_OccOrbToJ(occOrb) ...
                - obj.matpsi2.JK_OccOrbToK(occOrb);
            fockVec = reshape(obj.coreHamilt, [], 1) + reshape(gMat, [], 1);
        end
        
        function [orbital, orbEigValues] = SolveFockVec(~, fockVec, inv_S_Half)
            fockMat = reshape(fockVec, sqrt(numel(fockVec)), []);
            [orbitalOtho, orbEigValues] = eig(inv_S_Half*fockMat*inv_S_Half);
            [orbEigValues, ascend_order] = sort(diag(orbEigValues));
            orbital = inv_S_Half * orbitalOtho(:, ascend_order);
        end
        
        function elecEnergy = ElecEnergy(obj, fockVec, densVec)
            elecEnergy = (reshape(obj.coreHamilt, [], 1) + fockVec)'*densVec;
        end
        
        function newVec = Damping(~, dampingCoeff, vec, oldVec)
            coeffs = [dampingCoeff; (1 - dampingCoeff)];
            newVec = [vec, oldVec] * coeffs;
        end
        
        function cdiis = CDIIS(obj, numVectors)
            cdiis = CDIIS(obj.overlapMat, numVectors, 'r');
        end
        
        function ediis = EDIIS(obj, numVectors)
            ediis = EDIIS(obj.coreHamilt, numVectors, 'r');
        end
        
        function mciis = MCIIS(obj, numVectors)
            mciis = MCIIS(obj.overlapMat, numVectors, 'r');
        end
        
    end
    
    methods (Static)
                
        function properties = MatPsi2Interface(matpsi2)
            properties.overlapMat = matpsi2.Integrals_Overlap();
            properties.coreHamilt = matpsi2.Integrals_Kinetic() + matpsi2.Integrals_Potential();
            properties.nucRepEnergy = matpsi2.Molecule_NucRepEnergy();
            properties.numElectrons = matpsi2.Molecule_NumElectrons();
            properties.matpsi2 = matpsi2;
        end
        
    end
    
end