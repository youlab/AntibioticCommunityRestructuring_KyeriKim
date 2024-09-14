% Created by Kyeri Kim
% last edit: 4/4/2024 R2021a
close all
clear all
%% colors
red = [0.6350    0.0780    0.1840];
org = [0.8500    0.3250    0.0980];
ylw = [0.9290    0.6940    0.1250];
grn = [0.4660    0.6740    0.1880];
sky = [0.3010    0.7450    0.9330];
blu = [0         0.4470    0.7410];
ppl = [0.4940    0.1840    0.5560];

co1 = [247 150 70];
co11 = [247 222 197];
co2 = [75 172 198];
co22 = [201, 235, 245];

%% test LV formulation
% both negative and positive interaction, half and half
connectedness = 1;
neg_frac = 0.5;

%% random number set up
% rng shuffle
rng default

%% simulation parameter initializations
% basics
N 	= 30; % N subpopulations
x0 	= 1/100/N*ones(N,1); % mimic 1/100 dilution of initial growth phase, uniform initial density
tend 	= 8; % the length of simulation time, au.
% growth and lysis rates of subpopulations
max_Gs = linspace(0.7, 1.5, N)';
sigma = 1.50.*ones(N,1);
lambda = 0.1*ones(N,1);
max_inter = [0 max_Gs(end)*10 max_Gs(end)*20];

% G-L variations
L_vars = [0 0.1 0.2];

% pick random distribution numbers. These will be added to the sigma (slope) with the increased variations
m=1;s=2
while abs(m)>0.1 | s>1.2
random1 = randn(size(max_Gs),'like',max_Gs);
m=mean(random1);
s=std(random1);
end
% pick random distribution numbers. These will be added to the lambda (intercept) with the increased variations
m=1;s=2
while abs(m)>0.1 | s>1.2
random2 = randn(size(max_Gs),'like',max_Gs);
m=mean(random2);
s=std(random2);
end
sigmas = sigma+L_vars.*random1;
lambdas = lambda+L_vars.*random2;


%% Simulations and plots

% color setups for slow and fast subpopulations
global co sz alpha
alpha = .6;
coslow = [linspace(co1(1), co11(1), N/2)',linspace(co1(2), co11(2), N/2)',linspace(co1(3), co11(3), N/2)']/256;
cofast = [linspace(co22(1), co2(1), N/2)',linspace(co22(2), co2(2), N/2)',linspace(co22(3), co2(3), N/2)']/256;
co = [coslow;cofast];
set(groot, 'defaultAxesColorOrder',co);

% ------------------------------------------------------------------------------------------------------------------
figure(1);clf; % G-L plot with 0, 10, 20% variations
set(gcf,'units','normalized','position',[0.5,0.45,0.5*length(L_vars)/4,0.3])

tiledlayout(1,length(L_vars),'Padding','compact','TileSpacing','tight');
for i = 1:length(L_vars);
    nexttile(i)
        x = [0 5];
        for j = 1:length(x0)
        y = x.*sigmas(j,i)+lambdas(j,i);
        plot(x,y,'color',[co(j,:),0.5]); hold on
        end
        sz1 = (1-sigmas(:,i)).*max_Gs-lambdas(:,i);
        sz = 70*(normalize(sz1, 'range')+0.5); % normalize size for visualizations

    	scatter(max_Gs, max_Gs.*sigmas(:,i)+lambdas(:,i), sz, co, 'filled', 'MarkerEdgeColor',[0 0 0], 'MarkerFaceAlpha',alpha)
    	scatter(1.3,0.5,min(sz),'MarkerEdgeColor','k')
    	text(1.17,0.25,num2str(round(min(sz1),2)),'fontsize',8)
    	scatter(1.5,0.5,median(sz),'MarkerEdgeColor','k')
    	scatter(1.7,0.5,max(sz),'MarkerEdgeColor','k')
    	text(1.55,0.25,num2str(round(max(sz1),2)),'fontsize',8)        
   	axis square
        xticks([0:1:2])
        yticks([0:2:4])
        xticklabels({})
        yticklabels({})
        if i ==1
        xlabel('G')
        ylabel('L')%,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle');
        title('')
        xticklabels([0:1:2])
        yticklabels([0:2:4])
        else
        title([num2str(L_vars(i)*100),'%'])
        end
        ylim([0 4])
        xlim([0 2])

    ax = gca; 
    ax.FontSize = 15;
end


%%
figure(2);clf; % same G-L plot, but in G-L vs G format with 0, 10, 20% variations
set(gcf,'units','normalized','position',[0.5,0.45,0.5*length(L_vars)/4,0.3])

tiledlayout(1,length(L_vars),'Padding','compact','TileSpacing','tight');
for i = 1:length(L_vars);
    nexttile(i)
        x = [0 5];
        for j = 1:length(x0)
        y = (1-sigmas(j,i)).*x-lambdas(j,i);
        if i>1
        plot(x,y,'color',[co(j,:),0.3],'linewidth',1.1); hold on
        end
        plot(x, (1-1.5)*x-0.1, 'k-'); hold on
        end
        sz1 = (1-sigmas(:,i)).*max_Gs-lambdas(:,i);
        sz = 70*(normalize(sz1, 'range')+0.5);
        scatter(max_Gs, (1-sigmas(:,i)).*max_Gs-lambdas(:,i), sz, co, 'filled', 'MarkerEdgeColor',[0 0 0], 'MarkerFaceAlpha',alpha)

    	axis square
        xticks([0:1:2])
        yticks([-2:1:.2])
        xticklabels({})
        yticklabels({})
        
        if i ==1
        xlabel('G')
        ylabel('G - L');
        title('')
        xticklabels([0:1:2])
        yticklabels([-2:1:.2])
        end

        ylim([-2 0])
        xlim([0 2])

    ax = gca; 
    ax.FontSize = 15;
end

%% 3x3 plots of variations and interactions

c1 = 1/4*log(sum(x0)/sum(x0.*exp(4*max_Gs)));
c2 = 1/4*log(sum(x0.*exp(4*max_Gs))./sum(x0.*exp(4*max_Gs).*exp(4*((1-sigma).*max_Gs-lambda))));
xx = [-5 5];
yy = (1-sigma(1))*(xx-c1)-lambda(1)+c2;


figure(100);clf; % EF1 vs EF2, left to right: variation increase, top to bottom: max interaction increase

tiledlayout(length(max_inter),length(L_vars),'Padding','compact','TileSpacing','tight');
for i = 1:length(L_vars);
for j = 1:length(max_inter);
    s = (j-1)*length(L_vars)+i;
    ha(s) = nexttile(s);   
    plot(xx, yy, 'color', [1 1 1]*0.7, 'linewidth', 1.2); hold on
    xline(0, 'k--')
    yline(0, 'k--')
end
end

figure(101);clf; % Fraction vs Time, left to right: variation increase, top to bottom: max interaction increase
tiledlayout(length(max_inter),length(L_vars),'Padding','compact','TileSpacing','tight');
for i = 1:length(L_vars);
for j = 1:length(max_inter);
    s = (j-1)*length(L_vars)+i;
    hb(s) = nexttile(s);   
end
end


%% Model simulations

    % linear gLV model  
    max_neg = max_inter(end);
    max_pos = max_inter(end)*0.5; % scaled max positive interaction rate to avoid domination
    
    paramset = param_generator_lv(N,connectedness,neg_frac, max_neg, max_pos); % used code from Wu et al. (2022)
    
for j = 1:length(max_inter)
        params = paramset;
        params{1} = max_inter(j)/max_inter(end)*paramset{1};
    
for i = 1:length(L_vars)
    L_var = L_vars(i);
    
    s = (j-1)*length(L_vars)+i;
        sz = (1-sigmas(:,i)).*max_Gs-lambdas(:,i);
        sz = 70*(normalize(sz, 'range')+0.5);
   
    figure(100);
    [hf, hf2]= simplotfinal_linparam(x0,tend,params, max_Gs, sigmas(:,i), lambdas(:,i), 'off');
    hc  = get(hf,'children');
    hgc = get(hc, 'children');
    set(hgc, 'parent',ha(s));

    figure(101);
    hc2  = get(hf2,'children');
    hgc2 = get(hc2, 'children');
    set(hgc2, 'parent',hb(s));
end
end

% figure plot style editings
figure(101);
set(gcf,'units','normalized','position',[0,0.05,0.5*length(L_vars)/4,0.45/2*length(max_inter)])
for i = 1:length(L_vars);
for j = 1:length(max_inter);
    s = (j-1)*length(L_vars)+i;
    nexttile(s);
        xlim([0 8])
        ylim([0 1])
        xticks([0 4 8])
        yticks([0 0.5 1])
        xticklabels({})
        yticklabels({})
        axis square
        if s==length(L_vars)*(length(max_inter)-1)+1
        xlabel('Time (a.u.)')
        ylabel('Fraction')
        xticklabels([-4 0 4])
        yticklabels([0 0.5 1])
        end
    box off
    ax = gca; 
    ax.FontSize = 15;
    ax.XColor = 'black';
    ax.LineWidth =0.1;
end
end

%
figure(100)
set(gcf,'units','normalized','position',[0.5,0.05,0.5*length(L_vars)/4,0.45/2*length(max_inter)])
for i = 1:length(L_vars);
for j = 1:length(max_inter);
    s = (j-1)*length(L_vars)+i;
    nexttile(s);
        lims = 1.5;
        xlim([-lims lims]); ylim([-lims lims]);
        tix = [-lims -1 0 1 lims];
        tixx = [" ", -1, 0,1," "];
        xticks(tix);yticks(tix);
        xticklabels({})
        yticklabels({})
        axis square

        if s==length(L_vars)*(length(max_inter)-1)+1
        xlabel('EF_1')
        ylabel('EF_2')
        xticklabels(tixx);yticklabels(tixx);
        end
    box on
    ax = gca; 
    ax.FontSize = 15;
    ax.XColor = 'black';
end
end


%%===================================================================================================
function [fig1, fig2] = simplotfinal_linparam(x0,tend, params, max_Gs, sigma, lambda, figon)
global co sz alpha

[t,y] = run_timecourse_linparam(x0, tend, params, max_Gs, sigma, lambda);
fig1 = figure('visible', figon);

Nfrac = y./sum(y,2);
Ni = Nfrac(1,:);
t4 = find(t>=4, 1);
N0 = Nfrac(t4,:);
N4 = Nfrac(end,:);

N0Ni = log(N0./Ni)/4;
N4N0 = log(N4./N0)/4;

scatter(N0Ni, N4N0, sz, co,'filled', 'markeredgecolor',[0 0 0], 'MarkerFaceAlpha', alpha); hold on
axis square

fig2 = figure('visible', figon);
area(t,Nfrac, 'EdgeColor',[1 1 1]*0.95,'LineWidth',1.5,'FaceAlpha',0.9); hold on
xline(4, 'k--')
end

%-------------------------------------------------------------------------------------------------------
function param = param_generator_lv(N,connectedness,neg_frac,max_neg,max_pos,rnd_seed) % Wu et al. (2022)
% simplest function to generate the three sets of model parameters
% del, gam11, and gam22 all follow uniform distribution
% connectedness: the fraction of potential interactions that have
% interaction strengths (non-zero interaction strength)
% frac: the fraction of negative interaction over all interactions; This
% fraction is precise before zeroing out the diagnal of gam11 and gam22. 
% max_del: define the highest level of delta; minimum delta is 0.
% max_neg: define the maximum negative interaction strength; minimum
% negative interaction strength is 0
% max_pos: define the maximum positive interaction strength; minimum
% positive interaction strength is 0
% Set random number seed if a seed number is entered
if nargin==7
    rng(rnd_seed)
end

% Create a list of permutated indice for the interaction matrix
indsPerm = randperm(N*(N-1)); % exclude N diagonal values

% Initialize an interaction matrix
% An interaction is either positive or negative or neutral but cannot be
% both positive and negative. 
inter_init = zeros(N,N-1);

% Define the values of negative and positive interactions. 
neg_count = ceil(N*(N-1)*connectedness*neg_frac);
pos_count = floor(N*(N-1)*connectedness*(1-neg_frac));
gam = -rand(neg_count,1)*max_neg-0;
beta = +rand(pos_count,1)*max_pos+0;
%display(max_pos);

% Assign the negative and positive interaction strengths to the interaction
% matrix, using the permutated indice. 
inter_init(indsPerm(1:neg_count)) = gam;
inter_init(indsPerm((neg_count+1):(pos_count+neg_count))) = beta;

% add the 0 diaganol values to the interaction matrix
inter = zeros(N,N);
for i = 1:N
    if i==1
        inter(i,2:N) = inter_init(i,1:N-1);
    elseif i==N
        inter(i,1:N-1) = inter_init(i,1:end);
    else
        inter(i,1:i-1) = inter_init(i,1:i-1);
        inter(i,i+1:N) = inter_init(i,i:end);
    end
end

% Split the negative interaction matrix and positive interaction matrix
gam_f = inter.*(inter<0); % negative interactions
beta_f = inter.*(inter>0); % positive interactions

% Wrap together the three sets of parameters
param = {inter,gam_f,beta_f};
end

%-------------------------------------------------------------------------------------------------------
function [t,y] = run_timecourse_linparam(x0, tend, params, max_Gs, sigma, lambda);
    % start 4h growth without lysis
    [t1,y1] = run_core_ode_linparam(x0, 4, params, max_Gs, zeros(size(max_Gs)), zeros(size(max_Gs)));
    x0 = y1(end,:);
    % start growth and lysis
    [t2,y2] = run_core_ode_linparam(x0, tend-4, params, max_Gs, sigma, lambda);
    y = [y1;y2(2:end,:)];
    t = [t1;t2(2:end)+t1(end)];
end

%-------------------------------------------------------------------------------------------------------
function [t,y] = run_core_ode_linparam(x0, tend,params, max_Gs, sigma, lambda);
N = length(x0);
options=odeset('NonNegative',1:N,'AbsTol',1e-9); 
[t,y] = ode45(@core_ode_lv_linparam,[0 tend],x0,options,params{1},params{2},params{3}, max_Gs, sigma, lambda);
end
%-------------------------------------------------------------------------------------------------------
function dxdt = core_ode_lv_linparam(time, x, inter, gamma, beta, max_Gs,sigma, lambda);
dxdt = x.*((max_Gs.*(1-x) + inter*x).*(1-sigma) -lambda);
end