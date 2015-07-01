function energy = G09ReadEnergy(type)
if(nargin < 1)
    type = 'total';
end
if(strcmpi(type, 'NucRep'))
    keyword = ' nuclear repulsion energy ';
else
    keyword = ' E= ';
end
logFile = fopen('temp.log');
currLine = '';
while(ischar(currLine))
    if(~isempty(regexp(currLine, keyword, 'ONCE')))
        energy = str2double(regexp(currLine, '[+-]?[0-9]+.[0-9]+', 'match', 'ONCE'));
        break;
    end
    currLine = fgetl(logFile);
end
fclose(logFile);
end
