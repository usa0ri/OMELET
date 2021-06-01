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
    
%     f=f/max(f)*0.3; %normalize
    f=(f/max(f)).*0.05; %normalize
%     f=(f/max(f)).*0.03; %normalize
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
                fc(i,:),'FaceAlpha',alp,'EdgeColor',lc(i,:));
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
try
    set(gca,'TickLength',[0 0],'FontSize',8, 'Xtick',xtick_labels);
catch
    set(gca,'TickLength',[0 0],'FontSize',8, 'Xtick',[0.125 0.250]);
end
box on
if isempty(xL)==0
    set(gca,'XtickLabel',xL)
end

end