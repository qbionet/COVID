% Compute descriptive statistics of mobility, population, case, test data

clear all
close all
clc

%% Style choices
STYLE.colors.dark = {'#ef476f' '#ffd166' '#06d6a0' '#118ab2' '#adb5bd'};
STYLE.colors.light = {'#f7cad0' '#ffea00' '#b7e4c7' '#a2d2ff'};
STYLE.colors.anal  = {'#2d00f7' '#ff5400' '#f20089'};
thin  = 1;
thick = 5; 
fontSize = 12;

% time horizon, ticks and labels throughout the year
T_end = 360;
XLim_year = [0 366];
XTick_months = cumsum([0 31 29 31 30 31 30 31 31 30 31 30 31]);
XTickLabel_months = {'Jan' '' '' 'Apr' '' '' 'Jul' '' '' 'Oct' '' '' 'Jan'};
smoothWindow = 7;
miss_thr = Inf;  % maximum upscaling of case count

%% Load raw data
[parentdir,~,~]=fileparts(pwd);
RAW   = load([parentdir '/data_imported/data.mat']);
A         = RAW.A;
W         = RAW.W(:,:,1:T_end);
C         = smoothdata(RAW.C(:,1:T_end),2,'gaussian',smoothWindow);
T         = smoothdata(RAW.T(:,1:T_end),2,'gaussian',smoothWindow);
N         = RAW.N(:,1:T_end);
for i = 1:size(N,1)
    N(i,N(i,:)==0) = mean(N(i,:),2,'omitnan');
end
smoothdata(RAW.N(:,1:T_end),2,'gaussian',smoothWindow);
G         = RAW.G;
groups    = RAW.groups;
Z2G       = RAW.Z2G;
miss_rate = RAW.miss_rate;
ZIP_name  = RAW.ZIP_name;


%% Should small/heterophil ZIPs be combined?
COMBINE = 'Y';
if strcmp(COMBINE,'Y')
    targetName = 'ZIPcombined';
        
    % which ZIPs should be combined
    ZIP_list{1} = [53203 53233];
    ZIP_list{2} = [53205 53206]; 

    % composite ZIPs belong to these groups
    composite_group = [3 1]; 
    
    for i = 1:length(ZIP_list)
        ZIPs_to_combine = ZIP_list{i};
        [A, W, C, T, N, G, Z2G, ZIP_name] = ZIP_Combine(A, W, C, T, N, G, Z2G, ZIP_name, ZIPs_to_combine, composite_group(i));
    end         
else
    targetName = 'ZIPoriginal';
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

C(C<0) = 0;
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


%% Test upscaling
for t = 1:size(C,2)
    denom = 0;
    for g = 1:length(G)
        for i = 1:length(G{g})
            index = G{g}(i); 
            denom = denom + C(index,t)/group.testRate(g,t);
        end
    end
    D(t) = denom;
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

homophily_matrix = 100*sum(ZIP.W,3,'omitnan')./repmat(sum(ZIP.W,[2,3],'omitnan'),1,length(ZIP.name));

for i = 1:length(ZIP.name)
    ZIP.totalCaseRate(i) = 100*sum(ZIP.caseEst(i,:))/mean(ZIP.pop(i,:));
end



%% Fig S1 | poulation, movement, reported case rate
figure('Position',[0 0 1000 800])
tlt = tiledlayout(5, length(G));
tlt.Padding = 'none';
tlt.TileSpacing = 'compact';

% Population summary for each group
for g = 1:length(G)
    nexttile(), hold on
    for i = 1:length(G{g})
        plot(1e-4*ZIP.pop(G{g}(i),:),'LineWidth',thin,'Color',STYLE.colors.light{g})
    end
    grid on, xlim(XLim_year), ylim([0 15]), xtickangle(0)
    ylabel('population (10k)'), title(group.name{g})
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)  
end
for g = 1:length(G)
    nexttile()
    plot(1e-6*group.pop(g,:),'LineWidth',thick,'Color',STYLE.colors.dark{g})
    grid on, xlim(XLim_year), ylim([0 1]), xtickangle(0)
    ylabel('population (M)'), title(group.name{g})
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)   
end

% Trip rate summary
for g = 1:length(G)
    nexttile(), hold on
    for i = 1:length(G{g})
        plot(ZIP.tripRate(G{g}(i),:),'LineWidth',thin,'Color',STYLE.colors.light{g})
    end
    plot(group.tripRate(g,:),'LineWidth',thick,'Color',STYLE.colors.dark{g})
    grid on, xlim(XLim_year), ylim([0 20]), xtickangle(0)
    ylabel('trip rate'), title(group.name{g})
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months) 
end

% Reported cases for each group
for g = 1:length(G)
    nexttile(), hold on
    for i = 1:length(G{g})
        plot(100*ZIP.case(G{g}(i),:)./ZIP.pop(G{g}(i),:),'LineWidth',thin,'Color',STYLE.colors.light{g})
    end
    plot(100*group.case(g,:)./group.pop(g,:),'LineWidth',thick,'Color',STYLE.colors.dark{g})
    grid on, xlim(XLim_year), ylim([0 0.25]), xtickangle(0)
    ylabel('case rate (%)'), title(group.name{g})
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
end
for g = 1:length(G)
    nexttile(), hold on
    for i = 1:length(G{g})
        plot(100*cumsum(ZIP.case(G{g}(i),:))./ZIP.pop(G{g}(i),:),'LineWidth',thin,'Color',STYLE.colors.light{g})
    end
    plot(100*cumsum(group.case(g,:))./group.pop(g,:),'LineWidth',thick,'Color',STYLE.colors.dark{g})
    grid on, xlim(XLim_year), ylim([0 20]), xtickangle(0)
    ylabel('cumulative case rate (%)'), title(group.name{g})
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize) 
print(gcf,[parentdir '/figures/pop_mob_case_' targetName '.eps'],'-depsc')





%% Fig S2 | test and estimated case rate
figure('Position',[0 0 1000 800])
tlt = tiledlayout(5, length(G));
tlt.Padding = 'none';
tlt.TileSpacing = 'compact';

% Test data
for g = 1:length(G) 
   nexttile(), hold on
   plot(100*[squeeze(T(G{g},:))./N(G{g},:)]','LineWidth',thin,'Color',STYLE.colors.light{g})
   plot(100*group.testRate(g,:),'LineWidth',thick,'Color',STYLE.colors.dark{g})
   grid on, xlim(XLim_year), ylim([0 1.5]), xtickangle(0)
   ylabel('test rate (%)'), title(group.name{g})
   set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
end
for g = 1:length(G) 
   nexttile()
   plot(group.testScale(g,:),'LineWidth', thick,'Color',STYLE.colors.dark{g})
   grid on, xlim(XLim_year), ylim([1e0 1e3]), xtickangle(0)
   ylabel('case vs positive test'), title(group.name{g})
   set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months,'YScale','log')
end

% City vs county and testing
nexttile(), hold on
plot(1e-3*RAW.C_rep,'LineWidth',thick)
plot(1e-3*RAW.C_est,'LineWidth',thick)
grid on, xlim(XLim_year), ylim([0 15]), xtickangle(0)
ylabel('case count (1k)'), title('estimated and reported (county)')
set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)

nexttile()
plot(miss_rate,'LineWidth',thick)
grid on, xlim(XLim_year), ylim([1e0 1e3]), xtickangle(0)
ylabel('\xi'), title('estimated per reported (county)')
set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months,'YScale','log')

nexttile(), hold on
plot(1e-3*sum(C,1),'LineWidth',thick)
plot(1e-3*sum(C,1).*miss_rate(1:size(C,2))','LineWidth',thick)
grid on, xlim(XLim_year), ylim([0 20]), xtickangle(0)
ylabel('case count (1k)'), title('estimated and reported (city)')
set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)

nexttile(), hold on
plot(cumsum(1e-6*sum(C,1)),'LineWidth',thick)
plot(cumsum(1e-6*sum(C,1).*miss_rate(1:size(C,2))'),'LineWidth',thick)
grid on, xlim(XLim_year), ylim([0 1.2]), xtickangle(0)
ylabel('cumulative case count (1M)'), title('estimated and reported (city)')
set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)

% Case rate summary for each group
for g = 1:length(G)
    nexttile(3*length(G)+g), hold on
    for i = 1:length(G{g})
        plot(100*ZIP.caseEst(G{g}(i),:)./N(G{g}(i),:),'LineWidth',thin,'Color',STYLE.colors.light{g})
    end
    plot(100*group.caseEst(g,:)./group.pop(g,:),'LineWidth',thick,'Color',STYLE.colors.dark{g})
    grid on, xlim(XLim_year), ylim([0 1.5]), xtickangle(0)
    ylabel('case rate (%)'), title(group.name{g})
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)  
end
for g = 1:length(G)
    nexttile(4*length(G)+g), hold on
    for i = 1:length(G{g})
        plot(100*cumsum(ZIP.caseEst(G{g}(i),:)./N(G{g}(i),:),'omitnan'),'LineWidth',thin,'Color',STYLE.colors.light{g})
    end
    plot(100*cumsum(group.caseEst(g,:)./group.pop(g,:),'omitnan'),'LineWidth',thick,'Color',STYLE.colors.dark{g})
    grid on, xlim(XLim_year), ylim([0 80]), xtickangle(0)       
    ylabel('cumulative case rate (%)'), title(group.name{g})
    set(gca,'XTick',XTick_months,'XTickLabel',XTickLabel_months)
end
set(findall(gcf,'-property','FontSize'),'FontSize',fontSize)  
print(gcf,[parentdir '/figures/test_caseEst_' targetName '.eps'],'-depsc')



%% Save variables
close all
clearvars -except group ZIP city targetName comb ZIPs_to_combine parentdir
save([parentdir '/statistics/stats_' targetName '.mat'])











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         Function snippets                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Combine small ZIPs
function [A, W, C, T, N, G, Z2G, ZIP_name, comb, census] = ZIP_Combine(A, W, C, T, N, G, Z2G, ZIP_name, ZIPs_to_combine, composite_group)

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
    Z2G_comb(comb(1)) = composite_group;
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
end
