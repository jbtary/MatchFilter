% Create complete list of sac files/dates for all stations
% 
% INPUT:
% dirgen - path to directory with directories of all stations
% datemin - start date of dataset for complete list
% datemax - end date of dataset for complete list
% 
% OUTPUT:
% infoc - complete list of sac files
% datec - corresponding dates
% 
% Example:
% dirgen = '/media/jbtary/OSDisk/Users/jb.tary/Data/Marmara/2011/';
% datemin = datenum([2011,04,15,12,0,0]); % 15th of April, 12pm
% datemax = datenum([2011,07,31,23,0,0]); % 31st of July, 23am
% infoc = format_filelists(dirgen,datemin,datemax);

function [infoc,datec] = format_filelists(dirgen,datemin,datemax)

ndir = dir([dirgen 'G*']); ndir2 = dir([dirgen 'K*']);
ndir = [ndir;ndir2]; clear ndir2

% Complete list of sac files for each station
for ii = 1:size(ndir,1) % Loop on the number of stations
    
    info = dir([dirgen ndir(ii).name '/' ndir(ii).name 'Z/*.sac']);
    
    day = arrayfun(@(x) str2num(x.name(7:8)),info,'UniformOutput',false);
    month = arrayfun(@(x) str2num(x.name(10:11)),info,'UniformOutput',false);
    yr = arrayfun(@(x) str2num(x.name(13:16)),info,'UniformOutput',false);
    hr = arrayfun(@(x) str2num(x.name(18:19)),info,'UniformOutput',false);
    mm = arrayfun(@(x) str2num(x.name(21:22)),info,'UniformOutput',false);
    
    day = cell2mat(day); month = cell2mat(month); yr = cell2mat(yr);
    hr = cell2mat(hr); mm = cell2mat(mm);
    hrm = round(hr+(mm/60));
    datetmp = datenum([yr,month,day,hrm,zeros(size(hrm,1),1),zeros(size(hrm,1),1)]);
    
    [~,ia] = unique(datetmp);
    info = info(ia,:);
    datetmp = datetmp(ia,:);
    % To check which files were removed: setdiff(datetmp,datetmp2)
    
    eval(['date' ndir(ii).name '= datetmp;'])
    eval(['info' ndir(ii).name '= info;'])
    clear info day month yr hr mm hrm datetmp ia
end

% Exhaustive list of dates to compare with actual lists from OBSs

date(1,1) = datemin*24; datetmp = datemin*24; % To avoid roundoff errors
kk = 2;
while datetmp < datemax*24
    date(kk,1) = date(kk-1,1)+1;
    datetmp = date(kk,1);
    kk = kk+1;
end
date = date/24;

% Fit list of each OBS in complete file

infoc = zeros(size(date,1),size(ndir,1));
infoc = mat2cell(infoc,ones(1,size(date,1)),ones(1,size(ndir,1)));
datec = zeros(size(date,1),size(ndir,1));

for ii = 1:size(ndir,1)
    
    eval(['datetmp = date' ndir(ii).name ';'])
    eval(['info = info' ndir(ii).name ';'])
    
    [~,ia,ib] = intersect(datetmp,date); 
    name1 = extractfield(info,'name')';
    name1 = name1(ia,1);
    infoc(ib,ii) = name1;
    datec(ib,ii) = datetmp(ia,1);
    
    clear datetmp info ia ib name1
end
