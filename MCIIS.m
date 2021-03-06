classdef MCIIS < handle
    
    properties (Access = private)
        
        fockVectors;
        densVectors;
        S_Half;
        inv_S_Half;
        
        errorVectors;
        
    end
    
    methods
        
        function obj = MCIIS(overlapMatrix, numVectors, type)
            if(nargin < 2)
                numVectors = 5;
            end
            if(nargin < 3)
                type = 'r';
            end
            lenVec = numel(overlapMatrix);
            if(strcmpi(type, 'r'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.densVectors{1} = zeros(lenVec, numVectors);
                obj.errorVectors = zeros(lenVec, numVectors);
            elseif(strcmpi(type, 'u'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.densVectors{1} = zeros(lenVec, numVectors);
                obj.fockVectors{2} = zeros(lenVec, numVectors);
                obj.densVectors{2} = zeros(lenVec, numVectors);
                obj.errorVectors = zeros(2*lenVec, numVectors);
            end
            obj.S_Half = sqrtm(overlapMatrix);
            obj.inv_S_Half = obj.S_Half \ eye(size(obj.S_Half));
        end
        
        function Push(obj, newFockVector, newDensVector)
            for spin = 1:length(obj.fockVectors)
                % push Fock
                obj.fockVectors{spin}(:, 1:end-1) = obj.fockVectors{spin}(:, 2:end);
                obj.fockVectors{spin}(:, end) = newFockVector(:, spin);
                
                % push density
                obj.densVectors{spin}(:, 1:end-1) = obj.densVectors{spin}(:, 2:end);
                obj.densVectors{spin}(:, end) = newDensVector(:, spin);
            end
            
            % push new commutator error in
            obj.errorVectors(:, 1:end-1) = obj.errorVectors(:, 2:end);
            errorVec = [];
            for spin = 1:length(obj.fockVectors)
                FtDt = obj.inv_S_Half ...
                    * reshape(newFockVector(:, spin), sqrt(length(newFockVector(:, spin))), []) ...
                    * reshape(newDensVector(:, spin), sqrt(length(newDensVector(:, spin))), []) ...
                    * obj.S_Half;
                errorVec = [errorVec; reshape(FtDt - FtDt', [], 1)]; %#ok
            end
            obj.errorVectors(:, end) = errorVec;
        end
        
        function [optFockVector, coeffs, useFockVectors] = OptFockVector(obj)
            numVectors = sum(sum(obj.densVectors{1}.^2) ~= 0);
            [optFockVector, coeffs, useFockVectors] = obj.SolveForNumVectors(numVectors);
            while(isnan(coeffs))
                numVectors = numVectors - 1;
                [optFockVector, coeffs, useFockVectors] = obj.SolveForNumVectors(numVectors);
            end
        end
                
    end
    
    methods (Access = private)
        
        function [optFockVector, coeffs, useFockVectors] = SolveForNumVectors(obj, numVectors)
            useFockVectors = cell(1, length(obj.fockVectors));
            optFockVector = zeros(size(obj.fockVectors{1}, 1), length(obj.fockVectors));
            if(numVectors == 0 || numVectors == 1)
                for spin = 1:length(obj.fockVectors)
                    optFockVector(:, spin) = obj.fockVectors{spin}(:, end);
                    useFockVectors{spin} = obj.fockVectors{spin}(:, end);
                end
                if(length(useFockVectors) == 1)
                    useFockVectors = useFockVectors{1};
                end
                coeffs = 1;
                return;
            end
            useDensVectors = cell(1, length(obj.fockVectors));
            tensor = zeros(numVectors, numVectors, numVectors, numVectors);
            for spin = 1:length(obj.fockVectors)
                useFockVectors{spin} = obj.fockVectors{spin}(:, end-numVectors+1:end);
                useDensVectors{spin} = obj.densVectors{spin}(:, end-numVectors+1:end);
                tensor = tensor + obj.CalcTensor(useFockVectors{spin}, useDensVectors{spin});
            end
            
            useErrorVectors = obj.errorVectors(:, end-numVectors+1:end);
            onesVec = ones(numVectors, 1);
            hessian = [ ...
                useErrorVectors'*useErrorVectors, onesVec; ...
                onesVec', 0];
            diisCoefficients = hessian \ [zeros(numVectors,1); 1];
            iniCoeffs = diisCoefficients(1:end-1);
            
%             iniCoeffs = zeros(numVectors, 1);
%             iniCoeffs(end) = 1;
            
            coeffs = obj.NewtonSolver(tensor, iniCoeffs);
            
            for spin = 1:length(obj.fockVectors)
                optFockVector(:, spin) = useFockVectors{spin} * coeffs;
            end
            
            if(length(useFockVectors) == 1)
                useFockVectors = useFockVectors{1};
            end
%             disp(coeffs');
        end
        
        function tensor = CalcTensor(obj, useFockVectors, useDensVectors)
            nbf = sqrt(size(useFockVectors, 1));
            numVectors = size(useFockVectors, 2);
            error = zeros(numVectors^2, nbf^2);
            for dummyI = 1:numVectors
                fockI = reshape(useFockVectors(:, dummyI), nbf, []);
                for dummyJ = 1:numVectors
                    densJ = reshape(useDensVectors(:, dummyJ), nbf, []);
                    FtDt = obj.inv_S_Half * fockI * densJ * obj.S_Half;
                    error((dummyI-1)*numVectors+dummyJ, :) = reshape(FtDt - FtDt', [], 1);
                end
            end
            tensor = reshape(error*error', numVectors, numVectors, numVectors, numVectors);
        end
        
        function [value, grad, hess] = Target(~, coeffs, tensor)
            nVecs = length(coeffs);
            
            value = reshape(tensor, [], nVecs) * coeffs;
            value = reshape(value, [], nVecs) * coeffs;
            value = reshape(value, [], nVecs) * coeffs;
            value = reshape(value, [], nVecs) * coeffs;
            
            if(nargout > 1)
                grad = reshape( ...
                    permute(tensor, [1 2 3 4]) + ...
                    permute(tensor, [2 1 3 4]) + ...
                    permute(tensor, [3 1 2 4]) + ...
                    permute(tensor, [4 1 2 3]), [], nVecs) * coeffs;
                grad = reshape(grad, [], nVecs) * coeffs;
                grad = reshape(grad, [], nVecs) * coeffs;
            end
            
            if(nargout > 2)
                hess = reshape( ...
                    permute(tensor, [1 2 3 4]) + permute(tensor, [1 3 2 4]) + permute(tensor, [1 4 2 3]) + ...
                    permute(tensor, [2 1 3 4]) + permute(tensor, [2 3 1 4]) + permute(tensor, [2 4 1 3]) + ...
                    permute(tensor, [3 1 2 4]) + permute(tensor, [3 2 1 4]) + permute(tensor, [3 4 1 2]) + ...
                    permute(tensor, [4 1 2 3]) + permute(tensor, [4 2 1 3]) + permute(tensor, [4 3 1 2]), [], nVecs) * coeffs;
                hess = reshape(hess, [], nVecs) * coeffs;
                hess = reshape(hess, [], nVecs);
            end
            
        end
        
        function coeffs = NewtonSolver(obj, tensor, iniCoeffs)
            coeffs = iniCoeffs;
            for iter = 1:100
                [val, grad, hess] = obj.Target(coeffs, tensor);
                lambda = -4 * val;
                gradL = [grad + lambda; 0];
                onesVec = ones(length(coeffs), 1);
                hessL = [hess, onesVec; onesVec', 0];
%                 if(rcond(hessL) < 1e-20)
%                     disp('Inversion failed')
% %                     coeffs = iniCoeffs;
%                     return;
%                 end
                coeffsAndLambda = [coeffs; lambda];
                coeffsAndLambda = coeffsAndLambda - hessL \ gradL;
                coeffs = coeffsAndLambda(1:end-1);
                if(sqrt(mean(gradL.^2)) < 1e-12)
                    return;
                end
            end
            
            if(iter > 99)
                disp('Not converged');
%                 coeffs = iniCoeffs;
            end
        end
        
    end
    
end