classdef G09RSCF < RHF
    
    properties (Access = protected)
        
        info;
        
    end
    
    methods
        
        function obj = G09RSCF(info)
            G09RSCF.RunG09(info);
            matrices = G09RSCF.G09ReadMatrices({'overlap', 'coreHamilt'});
            properties.overlapMat = matrices{1};
            properties.coreHamilt = matrices{2};
            scalars = G09RSCF.G09ReadScalars({'NumElectrons', 'NucRepEnergy'});
            properties.numElectrons = scalars{1};
            properties.nucRepEnergy = scalars{2};
            properties.matpsi2 = [];
            obj@RHF(properties);
            obj.info = info;
        end
        
    end
    
    methods (Access = protected)
        
        function fockVec = OrbToFockVec(obj, orbital)
            obj.info.orbAlpha = orbital;
            G09RSCF.RunG09(obj.info);
            matrices = G09RSCF.G09ReadMatrices({'fockAlpha'});
            fockVec = reshape(matrices{1}, [], 1);
        end
        
        function elecEnergy = SCFEnergy(obj, fockVec, densVec)
            G09RSCF.RunG09(obj.info);
        end
        
        
        
    end
    
    methods (Static, Access = protected)
        
        scalars = G09ReadScalars(types);
        matrices = G09ReadMatrices(types);
        
        function RunG09(info)
            fileStr = G09RSCF.G09InputStr(info);
            gjfFile = fopen('temp.gjf', 'w');
            fprintf(gjfFile, '%s', fileStr);
            fclose(gjfFile);
            system('g09 temp.gjf');
            if(~G09RSCF.G09FileIsValid())
                throw(MException('G09RSCF:G09RSCF', 'g09 did not terminate correctly'));
            end
        end
        
    end
    
    methods (Static, Access = private)
        
        fileStr = G09InputStr(info);
        fileIsValid = G09FileIsValid();
        
    end
    
end
