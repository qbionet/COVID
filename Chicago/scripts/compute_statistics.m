% Compute descriptive statistics of mobility, population, case, test data

clear all
close all
clc

[parentdir,~,~]=fileparts(pwd);


%% Style choices
colors.dark = {'#ef476f' '#ffd166' '#06d6a0' '#118ab2' '#adb5bd'};
colors.light = {'#f7cad0' '#ffea00' '#b7e4c7' '#a2d2ff'};

colors.anal  = {'#2d00f7' '#ff5400' '#f20089'};
smoothWindow = 1;
miss_thr = Inf;  % maximum upscaling of case count
thin  = 1;
thick = 5; 
fontSize = 12;

% time horizon, ticks and labels throughout the year
T_end = 360;
XLim_year = [0 366];
XTick_months = cumsum([0 31 29 31 30 31 30 31 31 30 31 30 31]);
XTickLabel_months = {'Jan' '' '' 'Apr' '' '' 'Jul' '' '' 'Oct' '' '' 'Jan'};


%% Load raw data
RAW   = load([parentdir '/data_imported/data.mat']);
A         = RAW.A;
W         = RAW.W(:,:,1:T_end);
C         = RAW.C(:,1:T_end);
T         = RAW.T(:,1:T_end);
N         = RAW.N(:,1:T_end);
G         = RAW.G;
groups    = RAW.groups;
Z2G       = RAW.Z2G;
miss_rate = RAW.miss_rate;

ZIP_name  = RAW.ZIP_name;
census    = RAW.census; 


%% Should small/heterophil ZIPs be combined?
COMBINE = 'Y';
switch COMBINE
    case 'Y'
        ZIPs_to_combine = [60601 60602 60603 60604 60605 60606 60654 60661];
        [A, W, C, T, N, G, Z2G, ZIP_name, comb, census] = ZIP_Combine(A, W, C, T, N, G, Z2G, ZIP_name, ZIPs_to_combine, census);
        targetName = [parentdir '/figures/statistics/combined'];
        ZIP.census = census;
    case 'N'
        targetName = [parentdir '/figures/statistics/original'];
        ZIP.census = census;
end


%% Group names corrected
for g = 1:length(G)
    str = char(strrep(groups(g),'_',' '));
    str=lower(str);
    idx=regexp([' ' str],'(?<=\s+)\S','start')-1;
    str(idx)=upper(str(idx));
    group.name{g} = str;    
end

%% ZIP stats
for g = 1:length(G)
   for i = 1:length(G{g})
    ZIP.trip(G{g}(i),:) = squeeze(sum(W(G{g}(i),:,:),2));
   end
end
ZIP.pop      = N;
ZIP.tripRate = ZIP.trip./ZIP.pop;
ZIP.W = W;
ZIP.name = ZIP_name;
ZIP.area   = A;
ZIP.case = C;   % reported

%% Group stats
group.ZIPs = G;
for g = 1:length(G)
    % all [group] x [time]
    group.pop(g,:)   = sum(N(G{g},:),1);
    group.trip(g,:)  = sum(W(G{g},:,:),[1 2]);
    group.test(g,:)  = sum(T(G{g},:),1); 
    group.case(g,:)  = sum(C(G{g},:),1);
    group.area(g) = sum(ZIP.area(group.ZIPs{g}));

    group.tripRate(g,:) = group.trip(g,:)./group.pop(g,:);
    group.testRate(g,:) = group.test(g,:)./group.pop(g,:);
    group.caseRate(g,:) = group.case(g,:)./group.pop(g,:);
end

%% City stats
city.pop  = sum(group.pop,1);
city.trip = sum(group.trip,1);
city.tripRate = city.trip./city.pop;
city.Z2G  = Z2G;
city.area = sum(group.area);


%% City vs county and testing
figure('Position',[0 0 1000 250])
tlt = tiledlayout(1, 4);
tlt.Padding = 'none';
tlt.TileSpacing = 'compact';

nexttile
plot(1e-3*RAW.C_rep,'LineWidth',thick)
hold on
plot(1e-3*RAW.C_est,'LineWidth',thick)
axis square
grid on
xlim(XLim_year)
ylabel('case count (1k)')
set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
xtickangle(0)
title('estimated and reported (county)')

nexttile()
plot(miss_rate,'LineWidth',thick)
axis square
grid on
xlim(XLim_year)
ylim([1e0 1e3])
ylabel('\xi')
set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months,'YScale','log')
xtickangle(0)
title('estimated per reported (county)')

nexttile
plot(1e-3*sum(C,1),'LineWidth',thick)
hold on
plot(1e-3*sum(C,1).*miss_rate(1:size(C,2))','LineWidth',thick)
axis square
grid on
xlim(XLim_year)
ylabel('case count (1k)')
set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
xtickangle(0)
title('estimated and reported (city)')

nexttile
plot(cumsum(1e-6*sum(C,1)),'LineWidth',thick)
hold on
plot(cumsum(1e-6*sum(C,1).*miss_rate(1:size(C,2))'),'LineWidth',thick)
axis square
grid on
xlim(XLim_year)
ylabel('cumulative case count (1M)')
set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
xtickangle(0)
title('estimated and reported (city)')


set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)
if strcmp(COMBINE,'N')
    print(gcf,[targetName '/covidestim.eps'],'-depsc') 
end




%% Test upscaling
for t = 1:size(C,2)
    denom = 0;
    for g = 1:length(G)
        for i = 1:length(G{g})
            index = G{g}(i); 
            denom = denom + C(index,t)/group.testRate(g,t);
        end
    end
    case_estimate = miss_rate(t)*sum(C(:,t));
    xi0(t) = case_estimate/denom;
    
    for g = 1:length(G)
       group.testScale(g,t) = xi0(t)/group.testRate(g,t);
    end
end
group.caseEst = group.case.*group.testScale;
group.caseEst(isnan(group.caseEst)) = 0;

% ZIPs
for g = 1:length(G)
    for i = 1:length(G{g}) 
        ZIP.caseEst(G{g}(i),:) = group.testScale(g,:).*C(G{g}(i),:);
    end
end
ZIP.caseEst(isnan(ZIP.caseEst)) = 0;

% City
city.caseEst = sum(group.caseEst,1);


%% Census data
for g = 1:4
    index = group.ZIPs{g};
    group.census.pop(g)     = sum(ZIP.census.pop(index));
    group.census.hhold(g)   = ZIP.census.hhold(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.income(g)  = ZIP.census.income(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.hval(g)    = ZIP.census.hval(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.hyear(g)   = ZIP.census.hyear(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.hage(g)    = ZIP.census.hage(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.foreign(g) = ZIP.census.foreign(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.ptrans(g)  = ZIP.census.ptrans(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.health(g)  = ZIP.census.health(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.worker(g)  = ZIP.census.worker(index)'*ZIP.census.pop(index)/group.census.pop(g);   
    group.census.poverty(g) = ZIP.census.poverty(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.insur(g)   = ZIP.census.insur(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.ocrowd(g)  = ZIP.census.ocrowd(index)'*ZIP.census.pop(index)/group.census.pop(g);
    group.census.bo50(g)    = ZIP.census.bo50(index)'*ZIP.census.pop(index)/group.census.pop(g);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%       Movement and demographic plots         %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Population summary for each group
figure('Position',[0 0 1000 500])
tlt = tiledlayout(2, 4);
tlt.Padding = 'none';
tlt.TileSpacing = 'compact';

for g = 1:length(G)
    nexttile
    hold on
    for i = 1:length(G{g})
        plot(1e-4*ZIP.pop(G{g}(i),:),'LineWidth',thin,'Color',colors.light{g})
    end
    axis square
    grid on
    xlim(XLim_year)
    ylim([0 15])
    ylabel('population (10k)')
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
    xtickangle(0)
    title(group.name{g})
end
for g = 1:length(G)
    nexttile
    plot(1e-6*group.pop(g,:),'LineWidth',thick,'Color',colors.dark{g})
    axis square
    grid on
    xlim(XLim_year)
    ylim([0 1])
    ylabel('population (M)')
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
    xtickangle(0)
    title(group.name{g})
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)  
print(gcf,[targetName '/population.eps'],'-depsc') 


%% Trip rate summary
figure('Position',[0 0 1000 250])
tlt = tiledlayout(1, 4);
tlt.Padding = 'none';
tlt.TileSpacing = 'compact';

for g = 1:length(G)
    nexttile
    hold on
    for i = 1:length(G{g})
        plot(ZIP.tripRate(G{g}(i),:),'LineWidth',thin,'Color',colors.light{g})
    end
    plot(group.tripRate(g,:),'LineWidth',thick,'Color',colors.dark{g})
    axis square
    grid on
    xlim(XLim_year)
    ylim([0 20])
    ylabel('trip rate')
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
    xtickangle(0)
    title(group.name{g})
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize) 
print(gcf,[targetName '/tripRate.eps'],'-depsc')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%        Case and test rate plots         %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reported cases for each group
figure('Position',[0 0 1000 500])
tlt = tiledlayout(2, 4);
tlt.Padding = 'none';
tlt.TileSpacing = 'compact';

for g = 1:length(G)
    nexttile
    hold on
    for i = 1:length(G{g})
        plot(100*ZIP.case(G{g}(i),:)./ZIP.pop(G{g}(i),:),'LineWidth',thin,'Color',colors.light{g})
    end
    plot(100*group.case(g,:)./group.pop(g,:),'LineWidth',thick,'Color',colors.dark{g})
    axis square
    grid on
    xlim(XLim_year)
    ylim([0 0.25])
    ylabel('case rate (%)')
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
    xtickangle(0)
    title(group.name{g})
end
for g = 1:length(G)
    nexttile
    hold on
    for i = 1:length(G{g})
        plot(100*cumsum(ZIP.case(G{g}(i),:))./ZIP.pop(G{g}(i),:),'LineWidth',thin,'Color',colors.light{g})
    end
    plot(100*cumsum(group.case(g,:))./group.pop(g,:),'LineWidth',thick,'Color',colors.dark{g})
    axis square
    grid on
    xlim(XLim_year)
    ylim([0 20])
    ylabel('cumulative case rate (%)')
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
    xtickangle(0)
    title(group.name{g})
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize) 
print(gcf,[targetName '/caseRate_reported.eps'],'-depsc')



%% Test data
figure('Position',[0 0 1000 500])
tlt = tiledlayout(2, 4);
tlt.Padding = 'none';
tlt.TileSpacing = 'compact';

for g = 1:length(G) 
   nexttile
   plot(100*[squeeze(T(G{g},:))./N(G{g},:)]','LineWidth',thin,'Color',colors.light{g})
   hold on
   plot(100*group.testRate(g,:),'LineWidth',thick,'Color',colors.dark{g})
   grid on 
   axis square
   xlim(XLim_year)
   ylim([0 1.5])
   ylabel('test rate (%)')
   title(group.name{g})
   set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
   xtickangle(0)
end
for g = 1:length(G) 
   nexttile
   hold on
   plot(group.testScale(g,:),'LineWidth', thick,'Color',colors.dark{g})
   grid on 
   axis square
   xlim(XLim_year)
   ylim([1e0 1e3])
   ylabel('case vs positive test')
   title(group.name{g})
   set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months,'YScale','log')
   xtickangle(0)
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize) 
print(gcf,[targetName '/testRate.eps'],'-depsc')



%% Case rate summary for each group
figure('Position',[0 0 1000 500])
tlt = tiledlayout(2, 4);
tlt.Padding = 'none';
tlt.TileSpacing = 'compact';

for g = 1:length(G)
    nexttile
    hold on
    for i = 1:length(G{g})
        plot(100*ZIP.caseEst(G{g}(i),:)./N(G{g}(i),:),'LineWidth',thin,'Color',colors.light{g})
    end
    plot(100*group.caseEst(g,:)./group.pop(g,:),'LineWidth',thick,'Color',colors.dark{g})
    axis square
    grid on
    xlim(XLim_year)
    ylim([0 1.5])
    ylabel('case rate (%)')
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
    xtickangle(0)
    title(group.name{g})
end
for g = 1:length(G)
    nexttile
    hold on
    for i = 1:length(G{g})
        plot(100*cumsum(ZIP.caseEst(G{g}(i),:)./N(G{g}(i),:)),'LineWidth',thin,'Color',colors.light{g})
    end
    plot(100*cumsum(group.caseEst(g,:)./group.pop(g,:)),'LineWidth',thick,'Color',colors.dark{g})
    axis square
    grid on
    xlim(XLim_year)
    ylim([0 80])
    ylabel('cumulative case rate (%)')
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
    xtickangle(0)       
    title(group.name{g})
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)  
print(gcf,[targetName '/caseRate_estimated.eps'],'-depsc')  





%% Save variables
close all
clearvars -except group ZIP city COMBINE comb ZIPs_to_combine parentdir
switch COMBINE
    case 'Y'
        save([parentdir '/statistics/stats_combined.mat'])
    case 'N'
        save([parentdir '/statistics/stats_original.mat'])
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         Function snippets                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Combine small ZIPs
function [A, W, C, T, N, G, Z2G, ZIP_name, comb, census] = ZIP_Combine(A, W, C, T, N, G, Z2G, ZIP_name, ZIPs_to_combine, census)

    for i = 1:length(ZIPs_to_combine)
       comb(i) = find(ZIP_name==ZIPs_to_combine(i));     
    end
    size_before = length(ZIP_name);
    size_after  = size_before - length(comb) + 1;
    RIGHT = diag(ones(1,size_before));

    for i = 2:length(comb)
       RIGHT(comb(i),comb(1)) = 1; 
    end
    RIGHT(:,comb(2:end)) = [];
    LEFT = RIGHT';

    A_comb = LEFT*A;
    for t = 1:size(C,2)
        W_comb(:,:,t) = LEFT*W(:,:,t)*RIGHT;
        C_comb(:,t) = LEFT*C(:,t);
        T_comb(:,t) = LEFT*T(:,t);
        N_comb(:,t) = LEFT*N(:,t);
    end
    Z2G_comb = Z2G;
    Z2G_comb(comb(2:end)) = [];
    for g = 1:4
        G_comb{g} = find(Z2G_comb==g);
    end
    name_comb = ZIP_name;
    name_comb(comb(2:end)) = [];

    A         = A_comb;
    W         = W_comb;
    C         = C_comb;
    T         = T_comb;
    N         = N_comb;
    G         = G_comb;
    Z2G       = Z2G_comb;
    ZIP_name  = name_comb;
    
    census.name    = ZIP_name;
    AVE = diag(census.pop);
    for i = 2:length(comb)
       AVE(comb(1),comb(i)) = census.pop(comb(i)); 
    end
    AVE(comb(2:end),:) = [];
    census.pop     = sum(AVE,2);
    AVE = AVE./repmat(sum(AVE,2),1,size(AVE,2));
    
    
    census.hhold   = AVE*census.hhold;
    census.income  = AVE*census.income;
    census.hval    = AVE*census.hval;
    census.hyear   = AVE*census.hyear;
    census.hage    = AVE*census.hage;
    census.foreign = AVE*census.foreign;
    census.ptrans  = AVE*census.ptrans;
    census.health  = AVE*census.health;
    census.worker  = AVE*census.worker;
    census.poverty = AVE*census.poverty;
    census.insur   = AVE*census.insur;
    census.ocrowd  = AVE*census.ocrowd;
    census.bo50    = AVE*census.bo50;
    
    
end




