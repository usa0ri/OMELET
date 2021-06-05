function c = my_colors(col,m,is_othercolor)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% othercolor
if is_othercolor
    types = who('-file','colorData.mat');
    if nargin < 1
        c = types;
    end
    if isnumeric(col)
        col = char(types(col));
    end

    data = load('colorData.mat',col);
    if isempty(filenames(data))
        c = [];
    else
        c = interp1(linspace(0,1,size(data.(n),1)),data.(n),linspace(0,1,m),'cubic');
        c(c<0) = 0;
        c(c>1) = 1;
    end
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch col
    case 'blue'
        % '#0000FF'
        colvec = [0 0 1];
    case 'red'
        % '#FF0000'
        colvec = [1 0 0];
    case 'green'
        % '#00FF00'
        colvec = [0 1 0];
    case 'cyan'
        % '#00FFFF'        
        colvec = [0 1 1]; 
    case 'pink'
        % '#FF00FF'
        colvec = [1 0 1];
    case 'yellow'
        % '#FFFF00'
        colvec = [1 1 0];
    case 'black'
        % '#000000'
        colvec = [0 0 0];
    case 'white'
        % '#FFFFFF'
        colvec = [1 1 1];
        
    case 'matlab_blue'
        % '#0072BD'
        colvec = [0 0.4470 0.7410];
    case 'matlab_orange'
        % '#D95319'
        colvec = [0.8500 0.3250 0.0980];
    case 'matlab_yellow'
        % '#EDB120'
        colvec = [0.9290 0.6940 0.1250];
    case 'matlab_purple'
        % '#7E2F8E'
        colvec = [0.4940 0.1840 0.5560];
    case 'matlab_lblue'
        % '#4DBEEE'
        colvec = [0.3010 0.7450 0.9330];
    case 'matlab_red'
        % '#77AC30'
        colvec = [0.6350 0.0780 0.1840];      
end
c = repmat(colvec,m,1);

end

end