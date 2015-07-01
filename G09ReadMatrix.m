function matrix = G09ReadMatrix(matType)
beginning = cell(1, length(matType));
ending = cell(1, length(matType));
for iMat = 1:length(matType)
    if(strcmpi(matType{iMat}, 'overlap'))
        beginning{iMat} = ' *** Overlap *** ';
        ending{iMat} = ' *** Kinetic Energy *** ';
    elseif(strcmpi(matType{iMat}, 'coreHamilt'))
        beginning{iMat} = ' ****** Core Hamiltonian ****** ';
        ending{iMat} = ' SVDSVc ';
    elseif(strcmpi(matType{iMat}, 'fockAlpha'))
        beginning{iMat} = ' Fock matrix \(alpha\):';
        ending{iMat} = ' Fock matrix \(beta\):';
    elseif(strcmpi(matType{iMat}, 'fockBeta'))
        beginning{iMat} = ' Fock matrix \(beta\):';
        ending{iMat} = ' E= ';
    end
end

matrix = {};

logFile = fopen('temp.log');
currLine = '';
while(ischar(currLine))
    for iMat = 1:length(matType)
        blocks = {};
        iBlock = 0;
        if(~isempty(regexp(currLine, beginning{iMat}, 'ONCE')))
            readLine = fgetl(logFile);
            while(isempty(regexp(readLine, ending{iMat}, 'ONCE')))
                allMatches = regexp(readLine, '[+-]?[0-9]+.[0-9]+D[+-][0-9][0-9]', 'match');
                if(~isempty(allMatches)) % rowNum with doubles
                    numsInALine = zeros(1, numOfNums);
                    for iMatch = 1:length(allMatches)
                        numsInALine(iMatch) = str2num(allMatches{iMatch}); %#ok
                    end
                    rowNum = str2double(regexp(readLine, '[0-9]+', 'match', 'ONCE'));
                    currentBlock(rowNum, :) = numsInALine; %#ok
                    blocks{iBlock} = currentBlock; %#ok
                else % integers
                    intMatches = regexp(readLine, '[0-9]+', 'match');
                    numOfNums = length(intMatches);
                    iBlock = iBlock + 1;
                    currentBlock = zeros(0, numOfNums);
                end
                readLine = fgetl(logFile);
            end
            currMat = [blocks{:}];
            if(triu(currMat, 1) == 0)
                currMat = currMat + currMat' - diag(diag(currMat));
            end
            matrix{iMat} = currMat; %#ok
            
            currLine = readLine;
        end
        
    end
    currLine = fgetl(logFile);
end
fclose(logFile);

end
