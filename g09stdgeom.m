
numFiles = 7165;
molecules = cell(1, numFiles);
for fileInd = 1:numFiles
str = fileread(['./geoms/qm/qm-', sprintf('%04d', fileInd-1), '.out']);
mol = MoleculeOld(str);

newLine = sprintf('\n');
printCommand = sprintf('%s\n', '#p am1 scf(maxcycle=1) iop(5/13=1)');
printTitle = sprintf('%s\n', 'title');
printMol = sprintf('%d %d\n', [0 1]);
for iAtom = 1:size(mol.cartesian, 1)
    printMol = [printMol, sprintf('%d %35.25f %35.25f %35.25f\n', mol.cartesian(iAtom, :))]; %#ok
end

fileStr = [printCommand, newLine, printTitle, newLine, printMol, newLine, newLine];


gjfFile = fopen('temp.gjf', 'w');
fprintf(gjfFile, '%s', fileStr);
fclose(gjfFile);
system('g09 temp.gjf');



beginning = ' Standard orientation: ';
ending = ' ---------------------';

logFile = fopen('temp.log');
currLine = '';
numOfNums = 6;
while(ischar(currLine))
    blocks = {};
    iBlock = 1;
    if(~isempty(regexp(currLine, beginning, 'ONCE')))
        readLine = fgetl(logFile);
        readLine = fgetl(logFile);
        readLine = fgetl(logFile);
        readLine = fgetl(logFile);
        readLine = fgetl(logFile);
        while(isempty(regexp(readLine, ending, 'ONCE')))
            allMatches = regexp(readLine, '[+-]?[0-9]+(\.[0-9]+)?', 'match');
            if(~isempty(allMatches)) % rowNum with doubles
                numsInALine = zeros(1, numOfNums);
                for iMatch = 1:length(allMatches)
                    numsInALine(iMatch) = str2num(allMatches{iMatch}); %#ok
                end
                rowNum = str2double(regexp(readLine, '[0-9]+', 'match', 'ONCE'));
                currentBlock(rowNum, :) = numsInALine; %#ok
                blocks{iBlock} = currentBlock; %#ok
            end
            readLine = fgetl(logFile);
        end
        currMat = [blocks{:}];
        if(triu(currMat, 1) == 0)
            currMat = currMat + currMat' - diag(diag(currMat));
        end
        
        currLine = readLine;
    end
    currLine = fgetl(logFile);
end
fclose(logFile);

cart = currMat(:, [2, 4:end]);
molecules{fileInd} = Molecule(cart);
disp(fileInd);
end