%% (Activation Reactions)
% Calculating Validation Percent for Models with Added Reactions
%Part 1
% Define the Path of Model and Validation Excel Files
Modelfile_path ='C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\Toy_model.xlsx';
Validationfile_path ='C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\Toy_model_exp.xlsx';
% Validationfile_path ='C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\Toy_model_exp_semi.xlsx';% Semi-quantitative data

% Part 2
% Initail Parameters
Int_time = 40;
Steady_time = 40;
R_non_input_no = 3; % The first reaction No. after input reactions in model
% Part 3
% Derive the Source and Output Components for Adding Reactions
namepos = strfind(Modelfile_path,'.xls');
if rem (namepos,1)== 0
    namestr = Modelfile_path(1:namepos-1);
    namestr = cellstr(namestr);
else
    disp ('Error: please insert the path for model file with correct extension');
    return
end
[specID,~,reactionIDs,paramList,ODElist,~,~] = util.xls2Netflux(namestr,Modelfile_path);
source=specID;
[w,n,EC50,~,~,~] = paramList{:};
OUT_ACT = zeros (length(source), length(specID)); % OR Gate
% OUT_ACT = zeros (length(source), length(reactionIDs)); % AND Gate
%Part 4
% Determining temporary model file for adding reactions (copy of model file)
Modelfile_path_temp1 ='C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\Toy_model_temp.xlsx';
Control= Model_Output(Modelfile_path_temp1, Validationfile_path, Int_time, Steady_time);
% Control= Model_Output_semi(Modelfile_path_temp1, Validationfile_path, Int_time, Steady_time); % Semi-quantitative data
display (['Control =',num2str(Control)]); 

% Part 5
% Calculating model validation percent for each added reaction
for i=1:length(source)
    disp (['Act',num2str(i)]); 
    ss=source(i); 
    for j=1:length(specID) %OR Gate
        Add_reaction = strjoin([ss, '=>',specID(j) ]); %OR Gate Activation reaction
%                   Add_reaction = strjoin(['!',ss, '%=>%',specID(j) ]); %OR Gate Inhibitory reaction
%                   Add_reaction = strrep(Add_reaction, ' ', ''); %OR Gate Inhibitory reaction
%                   Add_reaction = strrep(Add_reaction, '%', ' '); %OR Gate Inhibitory reaction
        R = {'middle','r7',Add_reaction,0.5,1.4,0.5}; %OR Gate with default parameters
        xlswrite('C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\Toy_model_temp.xlsx',R ,2 ,'A9'); %OR Gate
%     for j=R_non_input_no:length(reactionIDs) % AND Gate
%                  Add_reaction = strjoin([ ss, '&',reactionIDs(j) ]); % AND Activation reaction
%                  Add_reaction = strjoin(['!',ss,'&',reactionIDs(j) ]); % AND Inhibitory reaction
%                  Add_reaction = strrep(Add_reaction, ' ', '');% AND Inhibitory reaction
%                  R = {'middle',['r',num2str(j)],Add_reaction,w(j),n(j),EC50(j)}; % AND Gate
%                  xlswrite('C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\Toy_model_temp.xlsx',R ,2 ,['A',num2str(j+2)]); % AND Gate
        OUT_ACT(i,j) = Model_Output(Modelfile_path_temp1, Validationfile_path, Int_time, Steady_time); 
%         OUT_ACT(i,j) = Model_Output_semi(Modelfile_path_temp1, Validationfile_path, Int_time, Steady_time); % Semi-quantitative data
        disp (['Act',num2str(i),'*',num2str(j),'=',num2str(OUT_ACT(i,j))]); 
        R= {nan,nan,nan,nan,nan,nan};%OR Gate
        xlswrite('C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\Toy_model_temp.xlsx',R ,2 ,'A9');%OR Gate
%                 Add_reaction = strjoin([reactionIDs(j) ]); % AND Gate
%                 R = {'middle',['r',num2str(j)],Add_reaction,w(j),n(j),EC50(j)}; % AND Gate
%                 xlswrite('C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\Toy_model_temp.xlsx',R ,2 ,['A',num2str(j+2)]); % AND Gate
    end
end
OR_Output= [['Control = ',num2str(Control)],specID'; specID,num2cell(OUT_ACT)];%OR Gate
xlswrite('C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\OR_gate.xlsx',OR_Output); %OR Gate
% AND_Output= [['Control = ',num2str(Control)],reactionIDs'; specID,num2cell(OUT_ACT)];%AND Gate
% AND_Output(:,2:R_non_input_no)=[]; %AND Gate
% xlswrite('C:\Users\ak2jj\Desktop\Structural Revision Method for Large-scale Signaling Networks\AND_gate.xlsx',AND_Output); %AND Gate
delete ('ODEfun.m');