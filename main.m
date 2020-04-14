function main

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The PDFSO problem is a new approach for fail-safe optimization of 
% structures, which takes into account the available information on the 
% probability of occurrence of different accidental scenarios. The results
% lead to a less conservative and more appropriate design compared to the
% traditional fail-safe optimization, as the actual data of each accidental
% situation are included in the optimization process.

% Input data must be modified in the function: defineInputParameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[d0,optPrm,modelPrm,limStatePrm,damConfPrm] = defineInputParameters;

auxiliaryParameters(modelPrm)

%--------------------------------------------------------------------------
options = optimoptions(@fmincon,'Algorithm','sqp','Display','iter','Diagnostics','on',...
            'DiffMinChange',optPrm.DiffMinChange,...
            'DiffMaxChange',optPrm.DiffMaxChange,...
            'TolCon',optPrm.TolCon,'TolFun',optPrm.TolFun,...
            'OutputFcn',@outfun);

[d,fval,exitflag,output]=fmincon(@(d)Objfun(d),d0,[],[],[],[],optPrm.lb,optPrm.ub,...
                                 @(d)confun(d,modelPrm,limStatePrm,damConfPrm),options)
%--------------------------------------------------------------------------
plotKey=1;
save('plotKey.mat','plotKey')
[c,ceq] = confun(d,modelPrm,limStatePrm,damConfPrm);
save('c.mat','c')

counter = 0;
nameFolder='opt_results_it';
writeOutputResults(nameFolder,exitflag,output)
writeConstraints(modelPrm,nameFolder,c)
writeResponses(modelPrm,nameFolder)
writeActiveConstr(modelPrm,optPrm,nameFolder,counter,fval,d,c)
writeOvercomeLimitStateAndConfig(damConfPrm,nameFolder)

end

function [d0,optPrm,modelPrm,limStatePrm,damConfPrm] = defineInputParameters

% initial value of design variables:
d0 = [1 1];

% number of design variables
nDV = length(d0);

% number of damaged configurations:
nDamages = 60;

% probability of occurrence of each damaged configuration
pDamages = [ones(1,20)'*0.025;ones(1,20)'*0.015;ones(1,20)'*0.010]; % must have nDamages components
% pDamages = ones(1,nDamages)'/nDamages; % all damaged configurations have the same probability

if abs(sum(pDamages)-1) > 1e-8 % two numbers equal within a tolerance.
   error('Error. The sum of probabilities must be equal to 1')
elseif length(pDamages) ~= nDamages
    error('Error. Length of pDamages must be equal to nDamages')
end

% target probability of failure:
pf = 0.05;

% number of Load Cases
nLC = 2;

% number of Limit-states:
nDcon = 2;

% In this example, Coefficients for damaged configurations are randomly generated
rng(1);

%--------------------------------------------------------------------------
% Add here the limit-states and Load Cases to be considered:
% h = Load Case
% k = Limit-state
% response <= responseMax    G safety region: responseMax - response >= 0     G:1 - response/responseMax
h=1;
k=1;
resp{h,1}{k,1} = @(d,Cf) (1/(Cf(1) + Cf(4)*d(1)*d(2) + Cf(5)*d(1)^2 + Cf(6)*d(2)^2));
respMax{h,1}{k,1} = 1/4;
limitStateFun{h,1}{k,1} = @(resp,respMax) 1 - resp/respMax;
Cfd0{h,1}{k,1} = ones(1,6);
CfDamages{h,1}{k,1} = lhsdesign(nDamages,6);

h=1;
k=2;
resp{h,1}{k,1} = @(d,Cf) (1/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(5)*d(1)^2 + Cf(6)*d(2)^2));
respMax{h,1}{k,1} = 1/4;
limitStateFun{h,1}{k,1} = @(resp,respMax) 1 - resp/respMax;
Cfd0{h,1}{k,1} = ones(1,6);
CfDamages{h,1}{k,1} = lhsdesign(nDamages,6);

h=2;
k=1;
resp{h,1}{k,1} = @(d,Cf) (1/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(4)*d(1)*d(2) + Cf(6)*d(2)^2));
respMax{h,1}{k,1} = 1/4;
limitStateFun{h,1}{k,1} = @(resp,respMax) 1 - resp/respMax;
Cfd0{h,1}{k,1} = ones(1,6);
CfDamages{h,1}{k,1} = CfDamages{1,1}{1,1};

h=2;
k=2;
resp{h,1}{k,1} = @(d,Cf) (1/(Cf(1) + Cf(2)*d(1) + Cf(3)*d(2) + Cf(4)*d(1)*d(2) + Cf(5)*d(1)^2));
respMax{h,1}{k,1} = 1/4;
limitStateFun{h,1}{k,1} = @(resp,respMax) 1 - resp/respMax;
Cfd0{h,1}{k,1} = ones(1,6);
CfDamages{h,1}{k,1} = CfDamages{1,1}{2,1};

%--------------------------------------------------------------------------
lb=[1e-03,1e-03];   % lower bound values
ub=[]; % upper bound values

TolCon=1e-3; % constraint tolerance
TolFun=1e-3; % objective function tolerance
DiffMinChange=0.001; % minimum change in design variables for finite difference gradients
DiffMaxChange=0.004; % maximum change in design variables for finite difference gradients
%--------------------------------------------------------------------------

optPrm = struct('lb',lb,'ub',ub,'TolCon',TolCon,'TolFun',TolFun,'DiffMinChange',DiffMinChange,'DiffMaxChange',DiffMaxChange);
modelPrm = struct('nDcon',nDcon,'nLC',nLC,'nDV',nDV);
limStatePrm = struct('resp',{resp},'respMax',{respMax},'limitStateFun',{limitStateFun},'Cfd0',{Cfd0},'CfDamages',{CfDamages});
damConfPrm = struct('nDamages',nDamages,'pDamages',pDamages,'pf',pf);

save('modelPrm.mat','modelPrm')
end

function auxiliaryParameters(modelPrm)

% plot the CDFs at each iteration 0=no 1=yes
plotKey=0;
save('plotKey.mat','plotKey')

% finite difference gradients counter
n_CONT = 0;  
save('n_CONT.mat','n_CONT')

% hitory of design variables and objective function
history_d_opt = [];
history_fval_opt = [];
save('history_d_opt.mat','history_d_opt');
save('history_fval_opt.mat','history_fval_opt');
%--------------------------------------------------------------------------
% History of Results
counter = 'iter';
nameFolder= 'opt_results_it';
writeTitleObjFunAndDesVar(modelPrm.nDV,nameFolder,counter)

Results_history.designVariables = [];
Results_history.c = [];
Results_history.respList = [];
Results_history.probList = [];
save('Results_history.mat','Results_history')

end

function fval=Objfun(d)

fval=d(1)+d(2);

end

function [c,ceq]=confun(d,modelPrm,limStatePrm,damConfPrm)

load('plotKey.mat')
extension='.sh'; 

load('n_CONT.mat')
n_CONT = n_CONT + 1;  
save('n_CONT.mat','n_CONT')

ceq = [];

%----------------------------
% INTACT MODEL
%----------------------------
[gd0] = evaluateLimitStateIntactModel(modelPrm,limStatePrm,d);

[c] = evaluateDeterministicConstraints(modelPrm,gd0);

%----------------------------
% DAMAGED CONFIGURATIONS
%----------------------------
[G]=evaluateLimitStateDamagedConf(modelPrm,limStatePrm,damConfPrm,d);

[x_CDFs,f_CDFs] = builtCDFs(modelPrm,damConfPrm,G);
save('x_CDFs.mat','x_CDFs')
save('f_CDFs.mat','f_CDFs')

[p_CDFs] = obtainProbabilityFromCDFs(modelPrm,x_CDFs,f_CDFs,plotKey);

[c] = evaluateProbabilisticConstraints(modelPrm,damConfPrm,c,p_CDFs);

%--------------------------------------------------------------------------
load('Results_history.mat')
load('response_current_it.mat')
load('probabilities_current_it.mat')
Results_history.designVariables = [Results_history.designVariables;d];
Results_history.c = [Results_history.c;c];
Results_history.respList = [Results_history.respList;response_current_it];
Results_history.probList = [Results_history.probList;probabilities_current_it];
save('Results_history.mat','Results_history')

end

function stop = outfun(d,optimValues,state)

    stop = false;

     switch state
         case 'init'
%              hold on
         case 'iter'
            load('modelPrm.mat')
            load('history_d_opt.mat');
            load('history_fval_opt.mat');
            history_fval_opt = [history_fval_opt; optimValues.fval];
            history_d_opt = [history_d_opt; d];
            save('history_d_opt.mat','history_d_opt');
            save('history_fval_opt.mat','history_fval_opt');
            i=length(history_fval_opt);
            
            nameFolder='opt_results_it';
            writeObjFunAndDesVar(modelPrm.nDV,nameFolder,history_d_opt(i,:),history_fval_opt(i),i)
  
         case 'done'
%              hold off
         otherwise
     end
end

%--------------------------------------------------------------------------

function [gd0] = evaluateLimitStateIntactModel(modelPrm,limStatePrm,d)

response_current_it = [];

for h=1:modelPrm.nLC
    for k=1:modelPrm.nDcon
        resp = limStatePrm.resp{h,1}{k,1}(d,limStatePrm.Cfd0{h,1}{k,1});
        response_current_it = [response_current_it,resp];
        gd0_h_k = limStatePrm.limitStateFun{h,1}{k,1}(resp,limStatePrm.respMax{h,1}{k,1});
        gd0{h,1}{k,1} = gd0_h_k;
    end
end

save('response_current_it.mat','response_current_it')
save('gd0.mat','gd0')

end

function [c] = evaluateDeterministicConstraints(modelPrm,gd0)

l=0;

for h=1:modelPrm.nLC
    for k=1:modelPrm.nDcon
        l=l+1;
        gd0_h_k = gd0{h,1}{k,1};
        % g>=0 --> c=-g
        c(l) = -gd0_h_k;
    end
end

end

function [G] = evaluateLimitStateDamagedConf(modelPrm,limStatePrm,damConfPrm,d)

for h=1:modelPrm.nLC
    for k=1:modelPrm.nDcon
        for g=1:damConfPrm.nDamages
            resp(g,:) = limStatePrm.resp{h,1}{k,1}(d,limStatePrm.CfDamages{h,1}{k,1}(g,:));
            G_h_k(g,:) = limStatePrm.limitStateFun{h,1}{k,1}(resp(g,:),limStatePrm.respMax{h,1}{k,1});
        end
        G{h,1}{k,1} = G_h_k;
        responseDamagedConf{h,1}{k,1} = resp;
    end
end

save('G.mat','G')
save('responseDamagedConf.mat','responseDamagedConf')

end

function [x_CDFs,f_CDFs] = builtCDFs(modelPrm,damConfPrm,G)

for h=1:modelPrm.nLC
    for k=1:modelPrm.nDcon
        G_h_k = G{h,1}{k,1};
        [x_prov,index] = sort(G_h_k,'ascend');

        f_prov = damConfPrm.pDamages(index);

        x_limitState_k = unique(x_prov); % sorted and non-repeated values of x_prov
        f_PDF_limitState_k = zeros(1,length(x_limitState_k))';
        for i=1:length(x_limitState_k)
             f_PDF_limitState_k(i) = sum(f_prov(x_limitState_k(i) == x_prov));
        end 

        f_CDF_limitState_k = cumsum(f_PDF_limitState_k);

        if length(x_limitState_k) == 1 % The CDF can not have only 1 value
            x_value=x_limitState_k;
            epsilon = 1e-8;
            x_limitState_k = [x_value-epsilon;x_value;x_value+epsilon];    
            f_CDF_limitState_k = [0;0.5;1];
        end     

        x_CDFs{h,1}{k,1} = x_limitState_k;
        f_CDFs{h,1}{k,1} = f_CDF_limitState_k;
    end
end

end

function [p_CDFs] = obtainProbabilityFromCDFs(modelPrm,x_CDFs,f_CDFs,plotKey)

% Entro en CDF y calculo el valor de probabilidad
xq = 0;

for h=1:modelPrm.nLC
    for k=1:modelPrm.nDcon
        x_limitState_k = x_CDFs{h,1}{k,1};
        f_CDF_limitState_k = f_CDFs{h,1}{k,1};
        [fq] = obtainInterpolatedValue(x_limitState_k,f_CDF_limitState_k,xq);
        p_CDFs{h,1}{k,1} = fq;
        
        if plotKey ==1
            plot_CDF(k,h,xq,fq,x_limitState_k,f_CDF_limitState_k)
        end
    end
end

end

function [fq] = obtainInterpolatedValue(x_limitState_k,f_CDF_limitState_k,xq)

fq = interp1(x_limitState_k,f_CDF_limitState_k,xq,'pchip'); % Previous neighbor interpolation
if xq < min(x_limitState_k)
    fq = 0 - min(x_limitState_k);
elseif xq >= max(x_limitState_k)
    fq = 1;
end

end

function [c] = evaluateProbabilisticConstraints(modelPrm,damConfPrm,c,p_CDFs)

l=length(c);
probabilities_current_it = [];

for h=1:modelPrm.nLC
    for k=1:modelPrm.nDcon
        l=l+1;
        p = p_CDFs{h,1}{k,1};
        probabilities_current_it = [probabilities_current_it,p];
        % p<=pf
        c(l) = p/damConfPrm.pf - 1;
        save('probabilities_current_it.mat','probabilities_current_it')
    end
end

end

%--------------------------------------------------------------------------

function plot_CDF(LS,LC,xq,fq,x_limitState_k,f_CDF_limitState_k)

scrsz = get(0,'ScreenSize');   
h=figure('Position',[1,1,scrsz(3)*6/6,scrsz(4)*3/3]);
FIGTemp=gcf;
FigName=strcat('CDF_LS',num2str(LS),'_LC',num2str(LC));

p_piecewise = stairs(x_limitState_k,f_CDF_limitState_k,'LineStyle','-','color','b','LineWidth',0.5);
hold on

for j_points=1:length(x_limitState_k)
    p1=plot(x_limitState_k(j_points),f_CDF_limitState_k(j_points),'o','MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r');
end

pinit = min(x_limitState_k);
pfin  = max(x_limitState_k);

pvector = linspace(pinit,pfin,1000);

[p_pchip] = plot_approximated_CDF(pvector,x_limitState_k,f_CDF_limitState_k);

% extend the axes limits
lb_factor = (1 - sign(pinit)*0.05);
ub_factor = (1 + sign(pfin)*0.05);

location_x = 0.15;
location_y = 0.70;

p2=plot([xq, xq],[0, 1],'LineStyle','--','color',[0.5,0.5,0.5],'LineWidth',0.3);
p3=plot([lb_factor*pinit, ub_factor*pfin],[fq, fq],'LineStyle','--','color',[0.5,0.5,0.5],'LineWidth',0.3);
p4=plot(xq,fq,'o','MarkerSize',2,'MarkerEdgeColor','g','MarkerFaceColor','g');

if fq<0 % truco para no poner probabilidad negativa en la leyenda
    fq=0;
end

[hleg1,hobj1]=legend([p_piecewise p1 p_pchip p4],'Piecewise CDF','Points from CDF','Interpolated CDF',strcat('$F_{\hat{G}_{j}}(0)=P[\hat{G}_{j} \le 0] =',num2str(fq),'$'),...
    'Location','east','Orientation','vertical');

textobj1 = findobj(hobj1, 'type', 'text');
set(textobj1, 'Interpreter', 'latex', 'fontsize', 10);
pos1=get(hleg1, 'Position');
set(hleg1, 'Position', [location_x location_y 0.30 0.18]);

axis([lb_factor*pinit, ub_factor*pfin,0,1])
xlabel('$\hat{G}_{j}$','FontSize',10,'Interpreter','latex')
ylabel('$F_{\hat{G}_{j}}$','FontSize',10,'Interpreter','latex')
set(gca,'FontSize',10)
ax = gca;
ax.XAxis.MinorTick = 'off';
ax.YAxis.MinorTick = 'off';
ax.XAxis.MinorTickValues = lb_factor*pinit:(ub_factor*pfin-lb_factor*pinit)/8:ub_factor*pfin;
ax.YAxis.MinorTickValues = 0:0.1:1;

% % yaxis horizontal
% ylh = get(gca,'ylabel');
% ylp = get(ylh, 'Position');
% set(ylh, 'Rotation',0, 'Position',ylp-0.04, 'VerticalAlignment','middle', 'HorizontalAlignment','right')

% Paper position (to save files)
a1 = 20;
a2 = 12.2;
PaperPosition = [0 0 a1 a2];

set(FIGTemp,'PaperUnits','centimeters','PaperPosition',PaperPosition)

FileFormat1 = 'epsc';
FileFormat2 = 'tiff';
% FileFormat3 = 'png';
print(FIGTemp,strcat('-d',FileFormat1),'-r300', FigName);
print(FIGTemp,strcat('-d',FileFormat2),'-r300', FigName);
% print(FIGTemp,strcat('-d',FileFormat3),'-r300', FigName);
% saveas(FIGTemp,FigName,'fig');

close(FIGTemp);

end

function [p_pchip] = plot_approximated_CDF(pvector,x_limitState_k,f_CDF_limitState_k)

fqvector=[];
for j_point = pvector

    xq = j_point;
    fq = interp1(x_limitState_k,f_CDF_limitState_k,xq,'pchip'); % Previous neighbor interpolation
      
    fqvector=[fqvector,fq];
    
%     p1=plot(j_point,fq,'o','MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','none');

end

p_pchip = plot(pvector,fqvector,'LineStyle','-','color',[0.5,0.5,0.5],'LineWidth',0.5);

end

%--------------------------------------------------------------------------

function writeTitleObjFunAndDesVar(nDV,nameFolder,counter)

fid=fopen(strcat(nameFolder,'.txt'),'a'); % write design variables and fval for each iteration
    fprintf(fid,'%9s %12s', ...
                counter,'Obj. Fun');
    for k=1:nDV
        fprintf(fid,' %10s', ...
                     strcat('d',num2str(k)));
    end
    fprintf(fid,'\n');
fclose(fid);

end

function writeObjFunAndDesVar(nDV,nameFolder,d,fval,n_CONT)

fid=fopen(strcat(nameFolder,'.txt'),'a');
    fprintf(fid,'%9d %12.4f', n_CONT, fval);
    for k=1:nDV
        fprintf(fid,' %10.4f', d(k));
    end 
    fprintf(fid,'\n');
fclose(fid);

end

function writeOutputResults(nameFolder,exitflag,output)

fid=fopen(strcat(nameFolder,'.txt'),'a');
    fprintf(fid,' \n')
fclose(fid);

fid=fopen(strcat(nameFolder,'.txt'),'a');
fprintf(fid,' \n %10s %16s \n', ...
                 'exitflag','constrviolation');
fprintf(fid,' %10i %16.4e \n', ...
             exitflag,output.constrviolation);             
fclose(fid);

end

function writeConstraints(modelPrm,nameFolder,c)

fid=fopen(strcat(nameFolder,'.txt'),'a');
    
    l=0;
    fprintf(fid,'\n');
    g=0;
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            fprintf(fid,'%36s ', strcat('c(',num2str(l),')'));
        end
    end
    %Damages
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            fprintf(fid,'%36s ', strcat('c(',num2str(l),')'));
        end
    end

    l=0;
    fprintf(fid,'\n');
    g=0;
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            fprintf(fid,'%36.4e ', c(l));
        end
    end
    %Damages
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            fprintf(fid,'%36.4e ', c(l));
        end
    end
    fprintf(fid,'\n');
    
fclose(fid);

end

function writeResponses(modelPrm,nameFolder)

load('n_CONT.mat')
load('Results_history.mat')

l=0;
fid=fopen(strcat(nameFolder,'.txt'),'a');
    fprintf(fid,' \n');
    g=0;
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            fprintf(fid,'%36s ', strcat('response_d',num2str(g),'_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')')); 
        end
    end
    % Damages
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            fprintf(fid,'%36s ', strcat('probability_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')')); 
        end 
    end
    fprintf(fid,'\n');
    ll=0;
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            ll=ll+1;
            response = Results_history.respList(n_CONT,ll);
            fprintf(fid,'%36.4e ', response);
        end
    end    
    ll=0;
    % Damages
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            ll=ll+1;
            prob = Results_history.probList(n_CONT,ll);
            fprintf(fid,'%36.4e ', prob);
        end
    end
    fprintf(fid,'\n\n');
fclose(fid);

end

function writeActiveConstr(modelPrm,optPrm,nameFolder,counter,fval,d,c)

load('n_CONT.mat')
load('Results_history.mat')

if counter ~= 0
    writeTitleObjFunAndDesVar(modelPrm.nDV,nameFolder,counter)
    writeObjFunAndDesVar(modelPrm.nDV,nameFolder,d,fval,n_CONT)
end
%--------------------------------------------------------------------------

matrixLimStateActViol = [];

% ACTIVE CONSTRAINTS
l=0;
fid=fopen(strcat(strcat(nameFolder,'.txt')),'a');
    fprintf(fid,'\n');
    
    fprintf(fid,'%30s \n','----- ACTIVE CONSTRAINTS -----');
    ll=0;
    g=0;
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            ll=ll+1;
            if abs(c(l)) <= 0.02
                response = Results_history.respList(n_CONT,ll);
                fprintf(fid,'%24s', strcat('constr_d0_LC',num2str(h),'_LS',num2str(k)));
                fprintf(fid,'%10s %12.4e %36s %12.4e', strcat('c(',num2str(l),')='),c(l),strcat('response_d',num2str(g),'_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')'),response);
                fprintf(fid,'\n');
            end
        end    
    end
    
    ll=0;
    % Damages
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            ll=ll+1;
            if abs(c(l)) <= 0.02
                matrixLimStateActViol = [matrixLimStateActViol; [h,k]];
                prob = Results_history.probList(n_CONT,ll);
                fprintf(fid,'%24s', strcat('constr_danos_LC',num2str(h),'_LS',num2str(k)));
                fprintf(fid,'%10s %12.4e %36s %12.4e', strcat('c(',num2str(l),')='),c(l),strcat('prob_failure_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')'),prob);
                fprintf(fid,'\n');
            end
        end    
    end   
    
    fprintf(fid,'\n\n');
fclose(fid);

% VIOLATED CONSTRAINTS
l=0;
fid=fopen(strcat(strcat(nameFolder,'.txt')),'a');
    fprintf(fid,'\n');
    fprintf(fid,'%30s \n','----- VIOLATED CONSTRAINTS -----');
    ll=0;
    g=0;
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            ll=ll+1;
            if c(l) > optPrm.TolCon
                response = Results_history.respList(n_CONT,ll);
                fprintf(fid,'%24s', strcat('constr_d0_LC',num2str(h),'_LS',num2str(k)));
                fprintf(fid,'%10s %12.4e %36s %12.4e', strcat('c(',num2str(l),')='),c(l),strcat('response_d',num2str(g),'_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')'),response);
                fprintf(fid,'\n');
            end
        end    
    end
    
    ll=0;
    % Damages
    for h=1:modelPrm.nLC
        for k=1:modelPrm.nDcon
            l=l+1;
            ll=ll+1;
            if c(l) > optPrm.TolCon
                matrixLimStateActViol = [matrixLimStateActViol; [h,k]];
                prob = Results_history.probList(n_CONT,ll);
                fprintf(fid,'%24s', strcat('constr_danos_LC',num2str(h),'_LS',num2str(k)));
                fprintf(fid,'%10s %12.4e %36s %12.4e', strcat('c(',num2str(l),')='),c(l),strcat('prob_failure_LC',num2str(h),'_LS',num2str(k),' (l=',num2str(l),')'),prob);
                fprintf(fid,'\n');
            end
        end    
    end   
    fprintf(fid,'\n\n');
fclose(fid);

save('matrixLimStateActViol.mat','matrixLimStateActViol')

end

function writeOvercomeLimitStateAndConfig(damConfPrm,nameFolder)

load('G.mat')
load('responseDamagedConf.mat')
load('matrixLimStateActViol.mat')

[row,~]=size(matrixLimStateActViol);

fid=fopen(strcat(nameFolder,'.txt'),'a');
    fprintf(fid,' \n%69s \n', '----- DAMAGED CONFIGURATIONS WHERE THE LIMIT-STATE IS OVERCOME ----- '); 
fclose(fid);

for i=1:row
    h = matrixLimStateActViol(i,1);
    k = matrixLimStateActViol(i,2);
    G_h_k = G{h,1}{k,1};
    resp = responseDamagedConf{h,1}{k,1};
    for g=1:damConfPrm.nDamages
        if G_h_k(g)<0
            G_viol = G_h_k(g);
            resp_viol = resp(g);
            ID_damageConfig_viol = g;
            p_viol = damConfPrm.pDamages(g);
            fid=fopen(strcat(nameFolder,'.txt'),'a');
                fprintf(fid,'\n %20s %5s %35s %25s %25s %25s ',strcat('LS=',num2str(k)),strcat('LC=',num2str(h)),strcat('constraintValue=',num2str(G_viol)),strcat('DamagedConfiguration=',num2str(ID_damageConfig_viol)),strcat('ResponseValue=',num2str(resp_viol)),strcat('ProbabilityValue=',num2str(p_viol))); 
            fclose(fid);
        end
    end
end

end

