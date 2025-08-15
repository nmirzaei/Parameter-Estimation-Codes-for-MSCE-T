clc
clear
warning('off', 'all')
coreID = getenv('coreID');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%User input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
user=input('Which cancer type? (breast, crc, thyroid)','s');
Sex = input('Which sex? (Female, Male)','s');
rndinpt= input('Do you want to reset the randomizer based on cpu ID (required for HPC)? (Yes=1, No=0)');
cohort = input('Please enter the cohort you want to use:');
estinpt = input('Is there a size file in the folder? (Yes=1, No=0)');
parallel = input('Use parallel? (true, false)');
iteration = input('Show iterations? (Yes=iter, No=off)','s');
if strcmp(user,'breast')
    Menarche_Age = input('Please enter the age at Menarche for the cohort');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Doubling rates and alpha3-beta3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(user,'crc')
    %the tumor proliferation rate for CRC
    bm = 2.572;
    N0 = 2E8; %Number of normal stem cells
elseif strcmp(user,'breast')
    %the tumor proliferation rate for BrC
    bm = 5.244;
    N0 = 1.74E10; % https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.aaf9011&file=tomasetti.sm.pdf
elseif strcmp(user,'thyroid')
    %Calculating the tumor proliferation rate for ThC
    bm = 1.458;
    N0 = 6.5E7; %https://www.science.org/action/downloadSupplement?doi=10.1126%2Fscience.1260825&file=tomasetti_sm_rev.pdf
else
    error('Wrong input for cancer type! Remember they are case sensitive.')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time step size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = input('Enter incidence age model fit restriction:');
t = 0:T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reading data and sorting by age
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strcat('data1_5Yr_Age_',user,'_',Sex,'.csv');
data = csvread(filename,1,0);
AgeSortData = sortrows(data,3);  %sort all rows based on the age column 
cohortmin = min(AgeSortData(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumSample = 200;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creating fully populated Cell_num_array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diagxsortData = sortrows(data,4);
diagxmin = min(data(:,4));
diagxmax = max(data(:,4));
j=1;
for i=1:size(AgeSortData,1)
    if diagxsortData(i,4)==diagxmin
        cell_num(j,:) = [diagxsortData(i,4) diagxsortData(i,5) diagxsortData(i,end-1) diagxsortData(i,end)];
        diagxmin = diagxmin+1;
        j=j+1;
    end
end
if estinpt==0
    years_to_estimate = cohortmin:diagxmax;
    S = size(cell_num,1);
    % Fit a model to the existing data
    mdl = fitlm(cell_num(:,1), cell_num(:,4),'interactions');
    
    % Perform extrapolation
    data_estimated = predict(mdl, years_to_estimate');
    
    % Data deviation from the regression line
    SD = sqrt(sum((data_estimated(end:-1:end-S+1,1)-cell_num(end:-1:1,4)).^2)/S);
    
    %producing data using normal distribution and the calculated deviation
    for i = 1:numel(data_estimated)
        data_estimated(i,:) = normrnd(data_estimated(i,:), SD);
    end
    
    diagxmin = min(data(:,4));
    diagxmax = max(data(:,4));
    j=1;
    
    %This loop is added to replace the size data we actually had with the
    %estimated ones. Since we don't want to estimate the sizes we already have
    %actual data for.
    for i=1:size(data_estimated,1)
        if years_to_estimate(i)>=cell_num(1,1)
            data_estimated(i) = cell_num(j,4);
            j=j+1;
        end
    end
    %Create the full data array
    Size_data = [(cohortmin:diagxmax)',(cohortmin+4:diagxmax+4)',data_estimated];
    T = array2table(Size_data);
    T.Properties.VariableNames(1:3) = {'year1','year2','size'};
    filename2 = strcat('size_',user,'.csv');
    writetable(T,filename2)
elseif estinpt==1
    filename2 = strcat('size_',user,'.csv');
    Size_data = csvread(filename2,1,0);
else
    error('wrong input!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cohortmin=cohort;
idx = find(AgeSortData(:,1)==cohortmin);
data1 = AgeSortData(idx,6);
Yeardx = AgeSortData(idx,4);
Age = AgeSortData(idx,3);
Agemin = min(Age);
idx1 = find(Size_data(:,1)==Yeardx(1,:));
idx2 = find(Size_data(:,1)==Yeardx(end,:));
Cell_num_array =[Size_data(idx1-Agemin:idx2,1) Size_data(idx1-Agemin:idx2,2) Size_data(idx1-Agemin:idx2,3)];    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lb_p = [0 0 0 0 0 0 0]; %lower bounds for the prior
ub_p = [1E-5 1E-2 1E-1 1 18 1 18]; %upper bounds for the prior
guess = prior_guesser(length(lb_p),lb_p',ub_p',NumSample,coreID,rndinpt);

options_fmincon = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...
    'ConstraintTolerance', 1e-6); % Tighten constraint tolerance for fmincon
% Define GA options, including the fmincon options
opts = optimoptions('ga', ...
    'HybridFcn', {@fmincon, options_fmincon}, ...
    'PopulationSize', NumSample, ...
    'InitialPopulationMatrix', guess', ...
    'UseParallel',parallel,...
    'Display','iter',...
    'MaxGenerations', 4E3, ...
    'UseParallel', false, ...
    'TolCon', 1e-6);
LB=[0 0 0 0 0 0 0]; %lower allowed parameter bounds
UB=[1E-5 1E-2 1E-1 1 18 1 18];  %upper allowed parameter bounds
%A and b are for linear constraints corresponding to alpha3 and
%beta3
A = [0 0 0 -1 0 1 0;0 0 0 0 -1 0 1];
b = [0;0];
[theta1,fval,exitflag,output,population,scores] = ga(@(p) ODECalc_hazard_SimpleBirth(p,T,data1,Cell_num_array(:,3),Age,N0,bm), 7, A,b,[],[],LB',UB',[],[],opts);
theta2 = [theta1 fval];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileName1 = sprintf('theta_output_%d.mat', str2num(coreID));
save(fileName1,'theta2')

