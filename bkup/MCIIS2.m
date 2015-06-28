classdef MCIIS2 < handle
    
    properties (Access = private)
        
        fockVectors;
        densVectors;
        S_Half;
        inv_S_Half;
        
    end
    
    methods
        
        function obj = MCIIS2(overlapMatrix, numVectors, type)
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
            elseif(strcmpi(type, 'u'))
                obj.fockVectors{1} = zeros(lenVec, numVectors);
                obj.fockVectors{2} = zeros(lenVec, numVectors);
                obj.densVectors{1} = zeros(lenVec, numVectors);
                obj.densVectors{2} = zeros(lenVec, numVectors);
            end
            obj.S_Half = sqrtm(overlapMatrix);
            obj.inv_S_Half = obj.S_Half \ eye(size(obj.S_Half));
        end
        
        function Push(obj, newFockVector, newDensVector)
            % push new Fock and density in
            for spin = 1:length(obj.fockVectors)
                obj.fockVectors{spin}(:, 1:end-1) = obj.fockVectors{spin}(:, 2:end);
                obj.fockVectors{spin}(:, end) = newFockVector(:, spin);
                
                obj.densVectors{spin}(:, 1:end-1) = obj.densVectors{spin}(:, 2:end);
                obj.densVectors{spin}(:, end) = newDensVector(:, spin);
            end
        end
        
        function [optFockVector, coeffs, useFockVectors] = OptFockVector(obj)
            useFockVectors = cell(1, length(obj.fockVectors));
            optFockVector = zeros(size(obj.fockVectors{1}, 1), length(obj.fockVectors));
            numVectors = sum(sum(obj.densVectors{1}.^2) ~= 0);
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
            tensorCell = cell(1, length(obj.fockVectors));
            for spin = 1:length(obj.fockVectors)
                useFockVectors{spin} = obj.fockVectors{spin}(:, end-numVectors+1:end);
                useDensVectors{spin} = obj.densVectors{spin}(:, end-numVectors+1:end);
                tensorCell{spin} = obj.CalcTensor(useFockVectors{spin}, useDensVectors{spin});
            end
            
            
            iniCoeffs = zeros(numVectors, 1);
            iniCoeffs(end) = 1;
            
            coeffs = iniCoeffs;
            for iter = 1:10
                [val, grad, hess] = obj.Target(coeffs, tensorCell);
                lambda = -4 * val;
                gradL = [grad + lambda; 0];
                if(norm(gradL) < 1e-10)
                    break;
                end
                onesVec = ones(length(coeffs), 1);
                hessL = [hess, onesVec; onesVec', 0];
                coeffsAndLambda = [coeffs; lambda];
                coeffsAndLambda = coeffsAndLambda - hessL \ gradL;
                coeffs = coeffsAndLambda(1:end-1);
            end
            
            for spin = 1:length(obj.fockVectors)
                optFockVector(:, spin) = useFockVectors{spin} * coeffs;
            end
            
            if(length(useFockVectors) == 1)
                useFockVectors = useFockVectors{1};
            end
%             disp(coeffs');
        end
        
    end
    
    methods (Access = private)
        
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
        
        function [value, grad, hess] = Target(~, coeffs, tensorCell)
            nVecs = length(coeffs);
            
            value = 0;
            for spin = 1:length(tensorCell)
                temp = reshape(tensorCell{spin}, [], nVecs) * coeffs;
                temp = reshape(temp, [], nVecs) * coeffs;
                temp = reshape(temp, [], nVecs) * coeffs;
                temp = reshape(temp, [], nVecs) * coeffs;
                value = value + temp;
            end
            
            if(nargout > 1)
                grad = zeros(nVecs, 1);
                for spin = 1:length(tensorCell)
                    tensor = tensorCell{spin};
                    temp = reshape( ...
                        permute(tensor, [1 2 3 4]) + ...
                        permute(tensor, [2 1 3 4]) + ...
                        permute(tensor, [3 1 2 4]) + ...
                        permute(tensor, [4 1 2 3]), [], nVecs) * coeffs;
                    temp = reshape(temp, [], nVecs) * coeffs;
                    temp = reshape(temp, [], nVecs) * coeffs;
                    grad = grad + temp;
                end
            end
            
            if(nargout > 2)
                hess = zeros(nVecs, nVecs);
                for spin = 1:length(tensorCell)
                    tensor = tensorCell{spin};
                    temp = reshape( ...
                        permute(tensor, [1 2 3 4]) + permute(tensor, [1 3 2 4]) + permute(tensor, [1 4 2 3]) + ...
                        permute(tensor, [2 1 3 4]) + permute(tensor, [2 3 1 4]) + permute(tensor, [2 4 1 3]) + ...
                        permute(tensor, [3 1 2 4]) + permute(tensor, [3 2 1 4]) + permute(tensor, [3 4 1 2]) + ...
                        permute(tensor, [4 1 2 3]) + permute(tensor, [4 2 1 3]) + permute(tensor, [4 3 1 2]), [], nVecs) * coeffs;
                    temp = reshape(temp, [], nVecs) * coeffs;
                    temp = reshape(temp, [], nVecs);
                    hess = hess + temp;
                end
            end
            
        end
        
    end
    
end