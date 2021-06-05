function plotElasticity(obj,savedir)

    mkdir(savedir);
    [enames, edata] = load_enames(obj);

    fnames = {'subpro','cofactor','allo'};
    for i=1:length(enames)
        plot_e(enames{i},edata{i},fnames{i},savedir);
    end

end


function [enames, edata] = load_enames(obj)

fluxnames = obj.model_data.X.rxn.rxn_names_include;
metnames = obj.model_data.X.met.met_names;
num_rc = obj.model_data.X.num.num_rc;
S = obj.model_data.out.S(:,obj.model_data.X.idx.idx_include);

a = obj.par.a;
b = obj.par.b;
e_cofactor = obj.par.e_cofactor;
e_allo = obj.par.e_allo;

fluxname_irrev = obj.model_data.X.rxn.rxn_names_irrev;
enames_subpro = {};
edata_subpro = [];
for i=1:num_rc
    if strcmp(fluxnames{i},'Cs')
        enames_subpro = [enames_subpro; ['AcetylCoA' '_' fluxnames{i} ];...
            ['OAA' '_' fluxnames{i} ]];
        edata_subpro = [edata_subpro, a(:,i), b(:,i)];
    elseif ismember(fluxnames{i},fluxname_irrev) || strcmp(fluxnames{i},'Glud1')
        enames_subpro = [enames_subpro; [metnames{S(:,i)<0} '_' fluxnames{i} ]];
        edata_subpro = [edata_subpro, a(:,i)];
    else
        enames_subpro = [enames_subpro; [metnames{S(:,i)<0} '_' fluxnames{i} ];...
            [metnames{S(:,i)>0} '_' fluxnames{i} ]];
        edata_subpro = [edata_subpro, a(:,i), b(:,i)];
    end
end

enames_cofactor = {'NAD+_Gpd1';'NADH_Gpd1';...
    'NAD+_Ldha';'NADH_Ldha';...
    'NAD+_Mdh2';'NADH_Mdh2';...
    'NAD+_Glud1';'NADH_Glud1';...
    'ADP_Pklr';...
    'ATP_Pcx';...
    'GTP_Pck1';...
    'FAD_Sdha'};
edata_cofactor = e_cofactor(:,[1:9 11:13]);
edata_cofactor(:,[2 4 6 8]) = -edata_cofactor(:,[2 4 6 8]);

enames_allo = {'AMP_i_Fbp1';'Citrate_a_Fbp1';...
    'ATP_i_Pklr';'F1,6P_a_Pklr';'Ala_i_Pklr';'Phe_i_Pklr';...
    'AcetylCoA_a_Pcx';...
    'ATP_i_Glud1';'ADP_a_Glud1';'GTP_i_Glud1';'Leu_i_Glud1'};
edata_allo = [e_allo(:,1:2) e_cofactor(:,10) e_allo(:,3:end)];
edata_allo(:,[1 3 5 6 8 10 11]) = -edata_allo(:,[1 3 5 6 8 10 11]);

enames = {enames_subpro,enames_cofactor,enames_allo};
edata = {edata_subpro,edata_cofactor,edata_allo};


end

function plot_e(enames_,edata_,fnames_,savedir_now)

col_list = repmat([0.6 0.6 0.6],length(enames_),1);

fig = figure('visible','off');
[h,L,MX,MED,bw] = my_violin(edata_,...
    'xlabel',enames_,...
    'facecolor',col_list);
ax = gca;
ax.XTickLabelRotation = 45;
ax.TickLabelInterpreter = 'none';
ax.YLabel.String = 'Elasticity coefficients';
ylim([min(-1,ax.YLim(1)) max(1,ax.YLim(2))]);

set(findobj(gca,'type','axes'),'FontSize',10,'FontName','SansSerif');

fig.PaperUnits = 'inches';
fig.PaperSize = [length(enames_)/3 4];
fig.PaperPosition = [0 0 length(enames_)/3 4];
% saveas( gcf, [ savedir_now '/elasticity_' fnames_ '.png' ] );
saveas( gcf, [ savedir_now '/elasticity_' fnames_ '.pdf' ] );
close all;


end

function [h,L,MX,MED,bw]=my_violin(Y,varargin)
%defaults:
%_____________________
xL=[];
fc=[1 0.5 0];
lc='k';
alp=0.5;
mc='k';
medc='r';
b=[]; %bandwidth
plotlegend=0;
plotmean=0;
plotmedian=1;
x = [];
%_____________________
%convert single columns to cells:
if iscell(Y)==0
    Y = num2cell(Y,1);
end

%get additional input parameters (varargin)
if isempty(find(strcmp(varargin,'xlabel')))==0
    xL = varargin{find(strcmp(varargin,'xlabel'))+1};
end
if isempty(find(strcmp(varargin,'facecolor')))==0
    fc = varargin{find(strcmp(varargin,'facecolor'))+1};
end
if isempty(find(strcmp(varargin,'edgecolor')))==0
    lc = varargin{find(strcmp(varargin,'edgecolor'))+1};
end
if isempty(find(strcmp(varargin,'facealpha')))==0
    alp = varargin{find(strcmp(varargin,'facealpha'))+1};
end
if isempty(find(strcmp(varargin,'mc')))==0
    if isempty(varargin{find(strcmp(varargin,'mc'))+1})==0
        mc = varargin{find(strcmp(varargin,'mc'))+1};
        plotmean = 1;
    else
        plotmean = 0;
    end
end
if isempty(find(strcmp(varargin,'medc')))==0
    if isempty(varargin{find(strcmp(varargin,'medc'))+1})==0
        medc = varargin{find(strcmp(varargin,'medc'))+1};
        plotmedian = 1;
    else
        plotmedian = 0;
    end
end
if isempty(find(strcmp(varargin,'bw')))==0
    b = varargin{find(strcmp(varargin,'bw'))+1}
    if length(b)==1
        disp(['same bandwidth bw = ',num2str(b),' used for all cols'])
        b=repmat(b,size(Y,2),1);
    elseif length(b)~=size(Y,2)
        warning('length(b)~=size(Y,2)')
        error('please provide only one bandwidth or an array of b with same length as columns in the data set')
    end
end
if isempty(find(strcmp(varargin,'plotlegend')))==0
    plotlegend = varargin{find(strcmp(varargin,'plotlegend'))+1};
end
if isempty(find(strcmp(varargin,'x')))==0
    x = varargin{find(strcmp(varargin,'x'))+1};
end
%
if size(fc,1)==1
    fc=repmat(fc,size(Y,2),1);
end


% Calculate the kernel density
i=1;
for i=1:size(Y,2)
    
    if all(Y{i}==0)
        F(:,i)=nan(100,1);
        U(:,i)=nan(100,1);
        MED(:,i)=nan;
        MX(:,i)=nan;
        bw(:,i)=nan;
        
    else
    
    if isempty(b)==0
        [f, u, bb]=ksdensity(Y{i},'bandwidth',b(i));
    elseif isempty(b)
        [f, u, bb]=ksdensity(Y{i});
    end
    
    f=(f/max(f)).*0.05; %normalize
    F(:,i)=f;
    U(:,i)=u;
    MED(:,i)=nanmedian(Y{i});
    MX(:,i)=nanmean(Y{i});
    bw(:,i)=bb;
    
    end
    
end

%Check x-value options
if isempty(x)
    x = zeros(size(Y,2));
    setX = 0;
else
    setX = 1;
    if isempty(xL)==0
        disp('_________________________________________________________________')
        warning('Function is not designed for x-axis specification with string label')
        warning('when providing x, xlabel can be set later anyway')
        error('please provide either x or xlabel. not both.')
    end
end

% Plot the violins
i=1;
for i=i:size(Y,2)
    if ~all(Y{i}==0)
    if isempty(lc) == 1
        if setX == 0
            h(i)=fill([F(:,i)+i;flipud(i-F(:,i))],[U(:,i);flipud(U(:,i))],...
                fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],...
                fc(i,:),'FaceAlpha',alp,'EdgeColor','none');
        end
    else
        if setX == 0
            h(i)=fill([F(:,i)+i/8;flipud(i/8-F(:,i))],[U(:,i);flipud(U(:,i))],...
                fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
        else
            h(i)=fill([F(:,i)+x(i);flipud(x(i)-F(:,i))],[U(:,i);flipud(U(:,i))],...
                fc(i,:),'FaceAlpha',alp,'EdgeColor',lc);
        end
    end
    hold on
    if setX == 0
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i)),...
                interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i)) ],...
                [MX(:,i) MX(:,i)],mc,'LineWidth',1);
        end
        if plotmedian == 1
            p(2)=plot([interp1(U(:,i),F(:,i)+i/8,MED(:,i)),...
                interp1(flipud(U(:,i)),flipud(i/8-F(:,i)),MED(:,i)) ],...
                [MED(:,i) MED(:,i)],medc,'LineWidth',1);
        end
    elseif setX == 1
        if plotmean == 1
            p(1)=plot([interp1(U(:,i),F(:,i)+i,MX(:,i))+x(i)-i,...
                interp1(flipud(U(:,i)),flipud(i-F(:,i)),MX(:,i))+x(i)-i],...
                [MX(:,i) MX(:,i)],mc,'LineWidth',1);
        end
        if plotmedian == 1
            p(2)=plot([interp1(U(:,i),F(:,i)+i,MED(:,i))+x(i)-i,...
                interp1(flipud(U(:,i)),flipud(i-F(:,i)),MED(:,i))+x(i)-i],...
                [MED(:,i) MED(:,i)],medc,'LineWidth',1);
        end
    end
    
    end
end

% Add legend if requested
if plotlegend==1 & plotmean==1 | plotlegend==1 & plotmedian==1
    
    if plotmean==1 & plotmedian==1
        L=legend([p(1) p(2)],'Mean','Median');
    elseif plotmean==0 & plotmedian==1
        L=legend([p(2)],'Median');
    elseif plotmean==1 & plotmedian==0
        L=legend([p(1)],'Mean');
    end
    
    set(L,'box','off','FontSize',14)
else
    L=[];
end

% Set axis
if setX == 0
%     axis([0.5 size(Y,2)+0.5, min(U(:)) max(U(:))]);
    axis([min(1/8-F(:))-1/16 max(F(:)+size(Y,2)/8) min(U(:)) max(U(:))]);
elseif setX == 1
    axis([min(x)-0.05*range(x) max(x)+0.05*range(x), min(U(:)) max(U(:))]);
end

for i=1:size(Y,2)
   xtick_labels(i) = mean([F(:,i)+i/8;flipud(i/8-F(:,i))] );
end
xtick_labels = linspace(xtick_labels(1),xtick_labels(end),size(Y,2));

% Set x-labels
% xL2={''};
% i=1;
% for i=1:size(xL,2)
%     xL2=[xL2,xL{i},{''}];
% end
set(gca,'TickLength',[0 0],'FontSize',8, 'Xtick',xtick_labels);
box on
if isempty(xL)==0
    set(gca,'XtickLabel',xL)
end

end
