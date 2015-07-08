classdef MoleculeOld < handle
    
    properties (SetAccess = private)
        
        cartesian;
        charge = 0;
        multiplicity = 1;
        
    end
    
    
    properties (Constant, Hidden = true)
        
        Bohr2Angstrom = 0.529177249;
        
    end
        
    methods
        
        function obj = MoleculeOld(arg1, charge, multiplicity)
            if(nargin > 1)
                obj.charge = charge;
            end
            if(nargin > 2)
                obj.multiplicity = multiplicity;
            end
            exception = MException('Molecule:Molecule','Input format wrong.');
            if(ischar(arg1)) % string
                arg1 = [arg1, char(10)];
                arg1 = obj.ReplaceDelimitersWithSemicolon(arg1);
                firstLine = regexp(arg1, '(.+?);', 'match', 'once');
                firstLine = obj.ReplaceAtomicSymbolsWithAtomicNumbers(firstLine);
                firstLine = str2num(firstLine(1:end-1)); %#ok
                if(length(firstLine) == 1 || length(firstLine) == 7) % zmatrix string
                    obj.cartesian = obj.ZMatrixToCartesian(obj.ZMatStrToZMatrix(arg1));
                elseif(length(firstLine) == 4) % cartesian string
                    obj.cartesian = obj.CartStrToCartesian(arg1);
                else
                    throw(exception);
                end
            elseif(ismatrix(arg1)) % matrix
                if(size(arg1,2) == 7) % zmatrix
                    obj.cartesian = obj.ZMatrixToCartesian(arg1);
                elseif(size(arg1,2) == 4) % cartesian
                    obj.cartesian = arg1;
                else
                    throw(exception);
                end
            else
                throw(exception);
            end
        end
        
        function molString = MoleculeString(obj)
            numAtoms = size(obj.cartesian,1);
            symbols = repmat(' ', numAtoms, 5);
            periodicTable = obj.PeriodicTable();
            for iatom = 1:numAtoms
                atomSymbol = periodicTable{obj.cartesian(iatom, 1)};
                symbols(iatom, 1:length(atomSymbol)) = atomSymbol;
            end
            formSpec = '%.10f';
            spaces = repmat('  ', numAtoms, 1);
            molString = [symbols, num2str(obj.cartesian(:, 2), formSpec), ...
                spaces, num2str(obj.cartesian(:, 3), formSpec), ...
                spaces, num2str(obj.cartesian(:, 4), formSpec), ...
                repmat(char(10), numAtoms, 1)];
            molString = reshape(molString', 1, []);
        end
        
        function numAtoms = NumAtoms(obj)
            numAtoms = size(obj.cartesian, 1);
        end
        
        function numElectrons = NumElectrons(obj)
            numElectrons = sum(obj.cartesian(:, 1)) - obj.charge;
        end
        
        function distMat = DistanceMatrix(obj)
            distMat = dist(obj.cartesian(:, 2:end)');
        end
        
        function Plot(obj)
            tempFileName = [tempname(), '.xyz'];
            tempFile = fopen(tempFileName, 'w');
            xyzStr = [num2str(size(obj.cartesian,1)), char(10), 'temp', char(10), obj.MoleculeString()];
            fprintf(tempFile, xyzStr);
            system(['jmol ', tempFileName]);
            system(['rm ', tempFileName]);
        end
        
    end
    
    methods (Access = private)
        
        function cartesian = ZMatrixToCartesian(~, zMatrix)
            xyz = zeros(size(zMatrix,1) + 2, 3);
            xyz(1:2, [3 1]) = [1 1; 0 1];
            for i = 2:size(zMatrix,1)
                ref = zMatrix(i, [2 4 6]);
                bond = zMatrix(i, 3);
                if(i == 2)
                    ref(2) = i - 2;
                    angle = 90;
                else
                    angle = zMatrix(i, 5);
                end
                if(i <= 3)
                    ref(3) = i - 3;
                    dihedral = 0;
                else
                    dihedral = zMatrix(i, 7);
                end
                xyz(i+2, :) = xyz(ref(1)+2, :);
                xyzRef21 = xyz(ref(2)+2, :) - xyz(ref(1)+2, :);
                xyzRef21 = xyzRef21 ./ norm(xyzRef21);
                xyzRef32 = xyz(ref(3)+2, :) - xyz(ref(2)+2, :);
                xyzRef32 = xyzRef32 - xyzRef21 * xyzRef32' .* xyzRef21;
                if(norm(xyzRef32) > 0)
                    xyzRef32 = xyzRef32 ./ norm(xyzRef32);
                end
                xyz(i+2, :) = xyz(i+2, :) + bond .* ( ...
                    cosd(angle) .* xyzRef21 + ...
                    sind(angle) .* (cosd(dihedral).*xyzRef32 ...
                    - sind(dihedral) .* cross(xyzRef21, xyzRef32)));
            end
            cartesian = [zMatrix(:, 1) xyz(3:end, :)];
        end
        
        function zmatrix = ZMatStrToZMatrix(obj, zMatStr)
            zMatStr = obj.ReplaceDelimitersWithSemicolon(zMatStr);
            indDelim = regexp(zMatStr, ';');
            if(length(indDelim)<2)
                zMatStr = [zMatStr(1:indDelim(1)-1), ' 0 0 0 0 0 0;', zMatStr(indDelim(1)+1:end)];
            elseif(length(indDelim)<3)
                zMatStr = [zMatStr(1:indDelim(1)-1), ' 0 0 0 0 0 0;', ...
                    zMatStr(indDelim(1)+1:indDelim(2)-1), ' 0 0 0 0;', zMatStr(indDelim(2)+1:end)];
            else
                zMatStr = [zMatStr(1:indDelim(1)-1), ' 0 0 0 0 0 0;', ...
                    zMatStr(indDelim(1)+1:indDelim(2)-1), ' 0 0 0 0;', ...
                    zMatStr(indDelim(2)+1:indDelim(3)-1), ' 0 0;', zMatStr(indDelim(3)+1:end)];
            end
            zMatStr = obj.ReplaceAtomicSymbolsWithAtomicNumbers(zMatStr);
            zmatrix = str2num(zMatStr); %#ok
        end
        
        function cartesian = CartStrToCartesian(obj, cartStr)
            cartStr = obj.ReplaceAtomicSymbolsWithAtomicNumbers(cartStr);
            cartesian = str2num(cartStr); %#ok
        end
        
        function string = ReplaceAtomicSymbolsWithAtomicNumbers(obj, string)
            periodicTable = obj.PeriodicTable();
            for number = 1:length(periodicTable)
                symbol = periodicTable{number};
                string = regexprep(string, [symbol,'[^[A-z]]'], [num2str(number),' ']);
            end
        end
        
        function string = ReplaceDelimitersWithSemicolon(~, string)
            allDelimiters = regexp(string, '[^A-z[\.]0-9 +-]','match');
            string = regexprep(string, ['[',strjoin(allDelimiters,''),']'], ';');
        end
        
        function periodicTable = PeriodicTable(~)
            periodicTable = { ...
                'H',                                                                                                 'He', ...
               'Li', 'Be',                                                             'B' , 'C' , 'N' , 'O' , 'F' , 'Ne', ...
               'Na', 'Mg',                                                             'Al', 'Si', 'P' , 'S' , 'Cl', 'Ar', ...
               'K' , 'Ca', 'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', ...
               'Rb', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe', ...
               'Cs', 'Ba',  ...
                           'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', ...
                                 'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', ...
               'Fr', 'Ra', ...
                           'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', ...
                                 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn','Uut', 'Fl','Uup', 'Lv','Uus','Uuo'};
        end
        
    end
    
end
