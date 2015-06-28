classdef ADIIS < handle
    
    properties (Access = private)
        
        fockVectors;
        densVectors;
        
    end
    
    methods
        
        function obj = ADIIS(initVector, numVectors)
            if(nargin < 2)
                numVectors = 5;
            end
            obj.fockVectors = zeros(numel(initVector), numVectors);
            obj.densVectors = zeros(numel(initVector), numVectors);
        end
        
        function Push(obj, newFockVector, newDensVector)
            % push in matrices or rows or columns
            newFockVector = reshape(newFockVector, [], 1);
            newDensVector = reshape(newDensVector, [], 1);
            
            % push new Fock in
            obj.fockVectors(:, 1:end-1) = obj.fockVectors(:, 2:end);
            obj.fockVectors(:, end) = newFockVector;
            
            % push new density in
            obj.densVectors(:, 1:end-1) = obj.densVectors(:, 2:end);
            obj.densVectors(:, end) = newDensVector;
        end
        
        function newFockVector = Interpolate(obj)
            [useFockVectors, coeffs] = obj.FockCoeffs;
            newFockVector = useFockVectors * coeffs;
        end
        
        function [useFockVectors, coeffs] = FockCoeffs(obj)
            numVectors = sum(sum(obj.densVectors.^2) ~= 0);
            if(numVectors == 0 || numVectors == 1)
                useFockVectors = obj.fockVectors(:, end);
                coeffs = 1;
                return;
            end
            useFockVectors = obj.fockVectors(:, end-numVectors+1:end);
            useDensVectors = obj.densVectors(:, end-numVectors+1:end);
            
            % compute errors wrt. the latest Fock or density
            errorFockVectors = useFockVectors - ...
                repmat(useFockVectors(:,end), 1, numVectors);
            errorDensVectors = useDensVectors - ...
                repmat(useDensVectors(:,end), 1, numVectors);
            
            % first order term and Hessian
            firstOrder = 2.*(useFockVectors(:, end)'*errorDensVectors)';
            hessian = errorFockVectors'*errorDensVectors;
            hessian = hessian + hessian'; % multiply Hessian by 2 and cancels numerical error
            
            options = optimoptions('quadprog', 'Algorithm', 'active-set', 'Display', 'off');
            coeffs = quadprog(hessian, firstOrder, ...
                -eye(numVectors), zeros(numVectors, 1), ...
                ones(1,numVectors), 1, ...
                [], [], ...
                [], options);
            
%             disp(coeffs');
        end
        
    end
    
    methods (Access = private)
        
        function num = NumVectors(obj)
            num = size(obj.fockVectors, 2);
        end
        
        function varFull = ReducedGradient(obj, hessian, firstOrder, iniPoint)
            indDep = 1;
            indAct = 1:obj.NumVectors();
            indAct = indAct(indAct~=indDep);
            varFull = iniPoint;
            constr = ones(1, obj.NumVectors());
            for iter = 1:2000
                varDep = varFull(indDep);
                varAct = varFull(indAct);
                constrDep = constr(indDep);
                constrAct = constr(indAct);
                
                grad = firstOrder + hessian*varFull;
                gradDep = grad(indDep);
                gradAct = grad(indAct);
                
                redGrad = gradAct - constrAct'/constrDep*gradDep;
                
                delVarAct = zeros(obj.NumVectors()-1,1);
                delVarAct(redGrad<0) = -redGrad(redGrad<0);
                delVarAct(varAct>0) = -redGrad(varAct>0);
                
                if(norm(delVarAct) < 1e-8)
                    break;
                end
                
                stepSize = 2 ./ (2 + iter);
                varActSim = varAct + stepSize .* delVarAct;
                delVarAct(varActSim<0) = 0;
                
                delVarDep = - constrDep \ constrAct * delVarAct;
                varDepSim = varDep + stepSize .* delVarDep;
                
                if(varDepSim < 0)
                    strPos = indAct(varActSim>0);
                    indDep = strPos(1);
                    indAct = 1:obj.NumVectors();
                    indAct = indAct(indAct~=indDep);
                    continue;
                end
                
                delVarFull = zeros(obj.NumVectors(), 1);
                delVarFull(indDep) = delVarDep;
                delVarFull(indAct) = delVarAct;
                
                varFull = varFull + stepSize .* delVarFull;
                
            end
        end
        
    end
    
end