classdef G09RSCF < RHF
    
    properties (Access = protected)
        
        info;
        
    end
    
    methods
        
        function obj = G09RSCF(info)
            fileStr = G09RSCF.G09InputStr(info.chargeMult, info.cartesian, info.method, info.basisSet);
            G09RSCF.RunG09(fileStr);
            matrices = G09RSCF.G09ReadMatrices({'overlap', 'coreHamilt'});
            properties.overlapMat = matrices{1};
            properties.coreHamilt = matrices{2};
            properties.numElectrons = G09RSCF.G09ReadNumElectrons();
            properties.nucRepEnergy = G09RSCF.G09ReadEnergy('NucRep');
            properties.matpsi2 = [];
            obj@RHF(properties);
            obj.info = info;
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            obj.info.orbAlpha = orbital;
            fileStr = G09RSCF.G09InputStr(obj.info);
            G09RSCF.RunG09(fileStr);
            matrices = G09RSCF.G09ReadMatrices({'fockAlpha'});
            fockVec = reshape(matrices{1}, [], 1);
        end
        
        function elecEnergy = SCFEnergy(obj, fockVec, densVec)
            fileStr = G09RSCF.G09InputStr(obj.info);
            G09RSCF.RunG09(fileStr);
        end
        
    end
    
    methods (Static, Access = protected)
        
        fileStr = G09InputStr(info);
        fileIsValid = G09FileIsValid();
        energy = G09ReadEnergy(type);
        matrices = G09ReadMatrices(types);
        numElectrons = G09ReadNumElectrons();
        
    end
    
    methods (Static, Access = private)
        
        function RunG09(fileStr)
            gjfFile = fopen('temp.gjf', 'w');
            fprintf(gjfFile, '%s', fileStr);
            fclose(gjfFile);
            system('g09 temp.gjf');
            if(~G09RSCF.G09FileIsValid())
                throw(MException('G09RSCF:G09RSCF', 'g09 did not terminate correctly'));
            end
        end
        
    end
    
end
