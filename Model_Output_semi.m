function [OUTPUT] = Model_Output_semi(Modelfile_path, Validationfile_path, Int_time, Steady_time)
% Calculate new model's validation percent as OUTPUT
% Part 1
% Delete the previously formed ODE to make sure it's rewritten
currentfolder= pwd;
if exist([currentfolder '\ODEfun.m'],'file') == 2
    delete('ODEfun.m');
end
%
% Part 2
% Parse out model's name (xls2Netflux needs it)
namepos = strfind(Modelfile_path,'.xls');
if rem (namepos,1)== 0
    namestr = Modelfile_path(1:namepos-1);
    namestr = cellstr(namestr);
else
    disp ('Error: please insert the path for model file with correct extension');
    return
end
%
% Part 3
% generate ODE from model spreadsheet
[specID,reactionIDs,~,paramList,ODElist,~,~] = util.xls2Netflux(namestr,Modelfile_path);
commandLine = util.exportODE3(specID,paramList,ODElist); %Model V1
% commandLine = util.exportODE2(specID,paramList,ODElist); %Model V2
util.textwrite('ODEfun.m',commandLine);
%
% Part 4
% set up simulation options
tspan_ss = [0 Steady_time]; % run out to ss
tspan_int = [0 Int_time]; % run out to ss
%
% Part 5
% Read the validation sheet & Identification of each main columns
[~, txt, ~] = xlsread(Validationfile_path);
IncodeColumn = cellfun(@(x)isequal(x,'Input Code'), txt(1,:))';
OutputColumn = cellfun(@(x)isequal(x,'Output'), txt(1,:))';
MeasurementColumn = cellfun(@(x)isequal(x,'Measurement'), txt(1,:))';
IDColumn =  cellfun(@(x)isequal(x,'ID'), txt(1,:))';
%
% Part 6
% remove rows without data
noData = cellfun(@(x)isequal(x,'No Data'), txt(1:end, MeasurementColumn));
txt(noData, :) = [];
noData = cellfun(@isempty, txt(1:end, MeasurementColumn));
txt(noData, :) = [];
%
% Part 7
% Having validation inputs in workspace for debugging & Extracting Data from Validation file
assignin('base', 'txt', txt);
validationIDs = txt(2:end, IDColumn);
inputCode = txt(2:end, IncodeColumn);
measurement = txt(2:end,MeasurementColumn);
outputSpec = txt(2:end, OutputColumn);
UninputCode = unique(inputCode);
%
% Part 8
% convert species and rxn names to integer values to map name and reaction
% ID's to numerical integers, allowing the input code to be evaluated
% directly
for k = 1:length(specID)
    if isempty(specID{k})
        disp (['specID ',num2str(k),' missing']);
    else
        eval([specID{k},' = ',num2str(k),';']);
    end
end
for i = 1:length(reactionIDs)
    if isempty(reactionIDs{i})
        disp (['reactionIDs ',num2str(i),' missing']);
    else
        eval([reactionIDs{i},' = ',num2str(i),';']);
    end
end
for j = 1:length(validationIDs)
    if isempty(validationIDs{j})
        disp (['validationIDs ',num2str(j),' missing']);
    else
        eval([validationIDs{j}, ' = ', num2str(j), ';']);
    end
end

%
% Part 9
% Set validation threshold change
thresh1 = 0.001;
thresh2 = 1e-3*thresh1;
% options = [];
options = odeset('RelTol',0.1*thresh1,'AbsTol',0.1*thresh2);
% threshold, Ryall et al., 2012 set to 0.001 for sensitivity analysis
inc1={'LH'};
inc2={'MH'};
inc3={'HH'};
dec1={'LL'};
dec2={'ML'};
dec3={'HL'};
noc={'NC'};
inc={'Inc'};
dec={'Dec'};
numMatching = 0; % number of predictions consistent with the qualitative literature species behavior
%
% Part 10
% Define the size of some variable changing in loops and find indices of output species
outputSpeciesIndex = zeros(1, length(measurement));
yStartL = cell(1, length(UninputCode));
yEndL = cell(1, length(UninputCode));
prediction = cell(1, length(inputCode));
predChange = cell(1, length(inputCode));
match = zeros(1, length(measurement));
c =zeros (length(UninputCode),1);
for k = 1:length(outputSpec)
    [~,outputSpeciesIndex(k)] = ismember(outputSpec{k},specID);
end
%
% Part 11
% loop over all validation simulations read from the excel sheet
for i = 1:length(UninputCode)
    %     disp(['Simulation # ', num2str(i), ' of ',num2str(length(UninputCode))]) % write the simulation number to the command line to track loop progress
    [w,n,EC50,tau,ymax,y0] = paramList{:}; % reset params
    % Initial (control) Simulation
    a= strfind (UninputCode{i},';');
    b= strfind (UninputCode{i},'w(');
    if length(a)> 1.1 && length(b)> 0.1
        eval([UninputCode{i}(1:a(1)), '%', UninputCode{i}(a(1)+1:end)]);
        c(i)=1;
    end
    rpar = [w;n;EC50];
    params = {rpar,tau,ymax,specID};
    [~,y] = ode15s(@ODEfun, tspan_int, y0, options, params);
    yStart = y(end,:)'; % use the "no input" steady state as control
    % evaluate validation conditions from excel sheet
    eval(UninputCode{i});
    % Main Simulation
    rpar = [w;n;EC50];
    params = {rpar,tau,ymax,specID};
    [~,y] = ode15s(@ODEfun, tspan_ss, y0, options, params);
    yEnd = y(end,:)';
    % Determine Change of Species' Activity after Stimulation
    yStartL{i} = yStart;
    yEndL{i} = yEnd;
end
% Part 12
% Determination of activity change in each experiment
for i=1:length(inputCode)
    idx = find(ismember(UninputCode, inputCode{i}));
    activityChange_fold = abs(real(yEndL{idx}(outputSpeciesIndex(i)))/real(yStartL{idx}(outputSpeciesIndex(i))));
    
    % Determine type of Changes
    if activityChange_fold > 1.01 && activityChange_fold < 2.05 % increase
        prediction{i} = 'LH';
        predChange{i} = num2str(activityChange_fold);
        if isequal(inc1,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1; %if the simulation matches the experimental validation put a 1 in the vector
        else
            match(i) = 0; %if the simulation does not match put a 0 in the matrix
        end
        if isequal(inc,measurement(i))
            match(i)=match(i)+1;
            numMatching = numMatching + 1;
        end
    elseif activityChange_fold >= 2.05 && activityChange_fold < 5.05
        prediction{i} = 'MH';
        predChange{i} = num2str(activityChange_fold);
        if isequal(inc2,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1; %if the simulation matches the experimental validation put a 1 in the vector
        else
            match(i) = 0; %if the simulation does not match put a 0 in the matrix
        end
        if isequal(inc,measurement(i))
            match(i)=match(i)+1;
            numMatching = numMatching + 1;
        end
    elseif activityChange_fold >= 5.05
        prediction{i} = 'HH';
        predChange{i} = num2str(activityChange_fold);
        if isequal(inc3,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1; %if the simulation matches the experimental validation put a 1 in the vector
        else
            match(i) = 0; %if the simulation does not match put a 0 in the matrix
        end
        if isequal(inc,measurement(i))
            match(i)=match(i)+1;
            numMatching = numMatching + 1;
        end
    elseif activityChange_fold < 0.99 && activityChange_fold > 0.495
        prediction{i} = 'LL';
        predChange{i} = num2str(activityChange_fold);
        if isequal(dec1,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1; %if the simulation matches the experimental validation put a 1 in the vector
        else
            match(i) = 0; %if the simulation does not match put a 0 in the matrix
        end
        if isequal(dec,measurement(i))
            match(i)=match(i)+1;
            numMatching = numMatching + 1;
        end
    elseif activityChange_fold <= 0.495 && activityChange_fold > 0.195
        prediction{i} = 'ML';
        predChange{i} = num2str(activityChange_fold);
        if isequal(dec2,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1; %if the simulation matches the experimental validation put a 1 in the vector
        else
            match(i) = 0; %if the simulation does not match put a 0 in the matrix
        end
        if isequal(dec,measurement(i))
            match(i)=match(i)+1;
            numMatching = numMatching + 1;
        end
    elseif activityChange_fold <= 0.195
        prediction{i} = 'HL';
        predChange{i} = num2str(activityChange_fold);
        if isequal(dec3,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1; %if the simulation matches the experimental validation put a 1 in the vector
        else
            match(i) = 0; %if the simulation does not match put a 0 in the matrix
        end
        if isequal(dec,measurement(i))
            match(i)=match(i)+1;
            numMatching = numMatching + 1;
        end
    else % no change
        prediction{i} = 'NC';
        predChange{i} = num2str(activityChange_fold);
        if isequal(noc,measurement(i))
            numMatching = numMatching + 1;
            match(i) = 1;
        else
            match(i) = 0;
        end
    end
    
end
%
% Part 13
% Generate the outputs
OUTPUT = numMatching/length(measurement)*100;