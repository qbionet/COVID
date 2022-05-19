% Import ZIP and city level data (movement, population, demographic, case, area) 
% W: mobility data
% N: population data
% C: case data
% T: test data
% Z2G: ZIP to group mapping
% G: group to ZIP mapping
% groups: group names
% C_rep: city level reported case count  
% C_est: city level estimated case count 
% miss_rate: daily ratio of the above two
clear all
close all
clc

[parentdir,~,~]=fileparts(pwd);

census = LoadCensus([parentdir '/data_raw/ZIP/census.xlsx']);
[A ZIP_name] = LoadArea([parentdir '/data_raw/ZIP/area.xlsx']);
[C_rep C_est miss_rate] = LoadCase_RepEst([parentdir '/data_raw/city/est_count.xlsx']);
[W N C T Z2G G groups] = LoadNodeData([parentdir '/data_raw/ZIP/']);
save([parentdir '/data_imported/data.mat']);



%% Load ZIP level data
function [W N C T Z2G G groups] = LoadNodeData(fileNameBase)
    
    % Movement data
    fileName = [fileNameBase 'movement.xlsx'];
    opts = detectImportOptions(fileName);
    %preview(fileName,opts);
    sheets = str2double(sheetnames(fileName));
    for i = 1:length(sheets)
        index = sheets(i);
        opts.Sheet = pad(num2str(index),3,'left','0');
        M = readtable(fileName,opts);
        W(:,:,index) = table2array(M);
    end

    % Population size data
    fileName = [fileNameBase 'popsize.xlsx'];
    opts = detectImportOptions(fileName);
    %preview(fileName,opts);
    sheets = str2double(sheetnames(fileName));
    for i = 1:length(sheets)
        index = sheets(i);
        opts.Sheet = pad(num2str(index),3,'left','0');
        M = readtable(fileName,opts);
        N(:,index) = table2array(M(:,2));
    end

    % Case count data
    fileName = [fileNameBase 'count.xlsx'];
    opts = detectImportOptions(fileName);
    %preview(fileName,opts);
    sheets = str2double(sheetnames(fileName));
    for i = 1:length(sheets)
        index = sheets(i);
        opts.Sheet = pad(num2str(index),3,'left','0');
        M = readtable(fileName,opts);
        C(:,index) = table2array(M(:,2));
        T(:,index) = table2array(M(:,4));
        
    end
        
    fileName = [fileNameBase 'zips_to_groups.xlsx'];
    opts = detectImportOptions(fileName);
    opts.DataRange = 'A1';
    M = readtable(fileName,opts);
    groups = table2cell(unique(M));
    for i = 1:length(groups)
        index = find(strcmp(string(groups(i,:)),table2cell(M)));
        Z2G(index,1) = i;
        G{i} = index;
    end
    
    % Replace NaNs
    C(isnan(C)) = 0;
    T(isnan(T)) = 0;
end

%% Load city estimate data
function [C_rep C_est miss_rate] = LoadCase_RepEst(fileName)
    opts = detectImportOptions(fileName);
    M = readtable(fileName,opts);
    
    C_rep = table2array(M(:,2));
    C_est = table2array(M(:,3));
    miss_rate = C_est./C_rep;
end

%% Load census data
function census = LoadCensus(fileName);   
    opts = detectImportOptions(fileName);
    M = readtable(fileName,opts);
    
    census.name    = str2double(table2array(M(:,1)));
    census.pop     = table2array(M(:,2));
    census.hhold   = table2array(M(:,4));
    census.income  = table2array(M(:,5));
    census.hval    = table2array(M(:,6));
    census.hyear   = table2array(M(:,7));
    census.hage    = table2array(M(:,8));
    census.foreign = table2array(M(:,9));
    census.ptrans  = table2array(M(:,10));
    census.health  = table2array(M(:,14));
    census.worker  = table2array(M(:,15));
    census.poverty = table2array(M(:,18));
    census.insur   = table2array(M(:,19));
    census.ocrowd  = table2array(M(:,20));
    census.bo50    = table2array(M(:,22));       
end

%% Load area data
function [A ZIP_name] = LoadArea(fileName)   
    opts = detectImportOptions(fileName);
    M = readtable(fileName,opts);
    
    ZIP_name = table2array(M(:,1));
    A = str2double(table2array(M(:,2)));
    
end
