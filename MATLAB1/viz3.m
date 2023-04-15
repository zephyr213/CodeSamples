function M = viz3(Temp,StateM,trange,dt,mypath,varargin)
% viz function to create an animation to show the temperature evolution of
% each element
% all input variable names are taken directly from the build.m script.
% Temp - temperature matrix
% StateM - overall size of sample
% trange - vecotr - time range; scalar - ending time step
% varargin - used as a vector to specify the plane of interest (currently
% only support planes perpendicular with x,y,z axes)

Temp = Temp - 273.15;

% time range of interest
if length(trange) > 1
    % input is a vector (time step range)
    tstart = trange(1);
    tend = trange(2);
else
    % use the all time steps
    tstart = 1;
    tend = trange;
end
tstepvec = tstart:tend;
tvec = dt*tstepvec;

% first determine elements in the plane of interest
if nargin == 8
    px = varargin{1};
    py = varargin{2};
    pz = varargin{3};
else
    px = 0;
    py = 0;
    pz = 1;
end

% find the elements in the plane of interest
[elist,elist_logi] = finde(mypath,px,py,pz);

minT = min(Temp(elist_logi,:),[],'all');
maxT = max(Temp(elist_logi,:),[],'all');

% colormap related
cmap = turbo(100);
npad = round(100/(maxT-minT)*minT);
cmap1 = [linspace(1,cmap(1,1),npad)',linspace(1,cmap(1,2),npad)',linspace(1,cmap(1,3),npad)'];
cmap = [cmap1;cmap];

% determine the size of matrix
allsize = size(StateM);
if px ~= 0
    %baseM = nan(1,allsize(2),allsize(3));
    baseM = nan(allsize(2),allsize(3));
    str_xlabel = 'z';
    str_ylabel = 'y';
    v = [-90,90];
elseif py ~=0
    %baseM = nan(allsize(1),1,allsize(3));
    baseM = nan(allsize(1),allsize(3));
    str_xlabel = 'z';
    str_ylabel = 'x';
    v = [-90,90];
elseif pz ~=0
    %baseM = nan(allsize(1),allsize(2),1);
    baseM = nan(allsize(1),allsize(2));
    str_xlabel = 'y';
    str_ylabel = 'x';
    v = [0,90];
end
clims = [0,max(Temp,[],'all')];
%clims = clims+273.15;
M(length(tvec)) = struct('cdata',[],'colormap',[]); % struct for movie
for i = 1:length(tvec)
    Tslice = Temp(tstepvec(i),elist_logi); % row vector
    for j = 1:numel(elist)
        if ~isnan(Tslice(j))
            x = mypath(elist(j),1);
            y = mypath(elist(j),2);
            z = mypath(elist(j),3);
            if px ~= 0
                %x = 1;
                baseM(y,z) = Tslice(j);
            elseif py ~=0
                %y = 1;
                baseM(x,z) = Tslice(j);
            elseif pz ~=0
                %z = 1;
                baseM(x,y) = Tslice(j);
            end
            %baseM(x,y,z) = Tslice(j);
        end
    end
    baseM(isnan(baseM)) = 0;

    tstring = sprintf('Time = %.2f s',tvec(i));
    if i == 1
        txt = annotation('textbox',[.65 .78 .3 .2],'String',tstring,'EdgeColor','none');
        myplot = imagesc(baseM,clims);
        myf = gcf;
        %myf.Visible = 'off';
    else
        txt.String = tstring;
        myplot.CData = baseM;
    end

    %imagesc(baseM,clims)
    colormap(cmap)
    c = colorbar;
    c.Ticks=[floor(minT):25:ceil(maxT)];
    c.Limits = [floor(minT),ceil(maxT) ];
    c.Label.String = "Temperature \circC";
    string_ti = "Nx=" + num2str(px) + "  Ny=" + num2str(py) + "  Nz=" + num2str(pz);
    title(string_ti)
    xlabel(str_xlabel)
    ylabel(str_ylabel)
    view(v)

    drawnow

    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3)+100, pos(4)+ti(2)+ti(4)];
    M(i) = getframe(ax,rect);

    %pause(0.001)


end
%myf.Visible = 'on';


end

function [elist,elist_logi] = finde(mypath,px,py,pz)
% this sub function is used to find the elements in the plane of interests

pv = [px,py,pz];  % assume the plane of interest is parallel to either x,y,z (only one non-zero value in px,py,pz)
tmp1 = mypath - pv; % mypath are all positive integers. So only those elements of interest will result in a 0 value in the result
elist_logi = any(tmp1 == 0, 2); % find the 0 value in each row
tmp3 = 1:length(mypath); % assume the rows are always larger than 3 (columns)
elist = tmp3(elist_logi);
end


