% build3.m - try to add assisted laser heating effect
% build4.m - try to add writing data to file to save memory

clear,clc
tic
% change the following parameters as needed
fprintf('Loading parameters ...')
% material parameters
Rho = 1125; % kg/m3
kx = 0.17; ky = 0.17; kz = 0.17; % W/mK, thermal conductivity along x,y,z
%C_p = 1467.5;   % J/kg*K, Specific Heat at 1 atm of Air

% processing parameters
h = 100; % W/m2K, coefficient of convection
c2k = 273.15;
Tb = 90 + c2k; % bed T
Tc = 40 + c2k; % chamber T
Te = 240 + c2k; % extrusion T
v = 50e-3; % extrusion speed 50 mm/s taken from Garon's notes.

% laser-related parameters
LaserSwitch = true; % switch variable to consider laser heating or not
laserS = 11*v;  % laser speed
laserR = 4e-4/2;    % m, laser focal radius, data from experiment d=403 micrometers
laserP = 1.5; % W, laser power
laserG = laserP/pi/laserR^2;    % laser irradiation W/m2      
ab = 0.2  ;   % laser abosorption
% set 0.2 for 10.60 micrometer IR laser as a test

if LaserSwitch
    fprintf('Laser on\n')
else
    fprintf('Laser off\n')
end


% simulation parameters
%Lx = 2e-3; %
%Ly = Lx/2;
Lx = 20e-3;%2e-3; unit meter
Ly = 0.4e-3;
Lz = 25.4e-3/10;%0.2e-3;   % 1 layer, assume layer height is 0.2mm
%Lz = Lz*2;
dt = 0.01; % timestep, s

% Init NB list

[mypath,nx,ny,nz,dx,dy,dz] = genpath(dt, Lx, Ly, Lz*2)  ;
StateM = zeros(nx,ny,nz);
nb_list = nan(nx*ny*nz,6);

% calculate some constants
% face vector [Ax,Ay,Az,Ax,Ay,Az]
Avec = [dy*dz, dx*dz, dx*dy];
Avec = [Avec, Avec];
%kvec = [kx,ky,kz,kx,ky,kz];
dvec = [dx,dy,dz,dx,dy,dz];
Bvec = Avec./dvec;
C1 = Rho * dx * dy * dz/dt;
C2 = kz*dx*dy/dz;

counter_element = 0;

nsteps = nx*ny*nz;
natoms = nx*ny*nz;

% compute some laser related properties:
% determine the laser movement speed in terms of the grids.
% assume laser head only moves in direction of either x or y.
lsvx = laserS/dx*dt; % laser speed on the grid in dt in x direction
lsvy = laserS/dy*dt; % laser speed on the grid in dt in y direction
lsdq = laserP*dt*ab;   % laser energy absorbed in a timestep
lsrx = laserR/dx;   % laser radius in # x grids
lsry = laserR/dy;   % laser radius in # y grids
lsx = 0;    % laser location
lsy = 0;    % laser location
dirx = 1;   %moving direction
diry = 0;

% ----- output data parameters
pr_interval = 1; % output data interval
fsteps = fopen('Temp.txt','wt');
fsteps2 = fopen('Temp.bin','w');
% -----


% assume 1 element deposited at each time step
%eNum = 0; % element number to be deposited. In current assumption, it's the same as istep. Keep it for future expansions.

%Temp = nan(nsteps, natoms);
Temp_pre = nan(1,natoms); % Temp of (N-1) step -- added for writing data to file 
Temp_new = nan(1,natoms); % Temp of N step -- added for writing data to file

% Initialization - set up initial condition with a single element deposited
eNum = 1; % current element count
cx = mypath(eNum,1);
cy = mypath(eNum,2);
cz = mypath(eNum,3);
StateM(cx,cy,cz) = 1;
nb_list(eNum,:) = 0; % init
nb_list(eNum,6) = -1;
tstep = 1;  % time step tracker
%Temp(tstep, eNum) = Te;  % array to log T evolution of all elements
Temp_pre(1,eNum) = Te;
% currently set to log at each time step, may need to change to record at
% every N step to decrease memory usage.
tloops = toc;

% -- write initial data
fprintf(fsteps,'%d ',tstep);
fprintf(fsteps,'%f ', Temp_pre);
fprintf(fsteps, '\n');

intfmt = 'uint16';
fwrite(fsteps2,natoms,intfmt); % write the total number of element to the binary file
fwrite(fsteps2,tstep,intfmt);
fwrite(fsteps,Temp_pre,'single');
% -----
fprintf('Starting Simulation ...')
for istep = 2:natoms
    if mod(istep,50) == 0 % print simulation progress
        fprintf('Now: %d/%d ~ %.2f%%\n',istep,natoms,istep/natoms*100)
    end
    %T_now = Temp(tstep,:); % slice of current T
    T_now = Temp_pre;
    
    % determine the laser irradiated elements in one step
    total_Qi = zeros(1,natoms);
    % first identify all "top elements"
    elements_ind = 1:height(nb_list);
    topelements = (nb_list(:,3) == 0); % identify "top elements" neighbour3 == 0 ;
    topelements_xyz = mypath(topelements,1:2);
    topelements_a = zeros(height(topelements_xyz),1);
    nsteps_ls = 11;
    for ilaser = 1:nsteps_ls
        % affected elements
        % d = sqrt(sum((topelements_xyz - [lsx,lsy]).^2,2));
        d = abs(topelements_xyz(:,1)-0.5-lsx); % 1d distance for now
        topelements_a = topelements_a + (d<=(lsrx+1e-5)); % number of times covered by laser during its movement within 1 timestep
        % advance laser location
        lsx = lsx + dirx*1;
        if (lsx > nx) || (lsx < 0)
            dirx = -dirx;
            lsx = lsx + dirx*2;
        end
    end
    topelements_ratio = topelements_a/nsteps_ls; % percentage of irradiation in a time step.
    topelements_Qi = laserG*dx*dy*ab*topelements_ratio;
    total_Qi(topelements) = topelements_Qi;
    % end laser processing


    % first evolv temperature of exisiting elements
    Tnew = nan(1,eNum);
    parfor ia = 1:eNum  % use parallel for-loop to go over all existing elements
        % it seems faster to not use parfor on apple m1 (Possibly for small
        % systems only)
        mynblist = nb_list(ia,:);
        cc = mynblist > 1; % conduction contacts
        vc = mynblist == 0; % convection contacts

        %   Tnow = Temp(tstep,ia);
        Tnow = T_now(ia);   % T of current element

        [know,Cnow] = fitkc(Tnow);   % DEBUG: the interp1 function takes most of the computation time.
        if mynblist(6) == -1 % with bed conduction
            Qcb = C2*(Tb - Tnow);   % Heat conduction from bed
        else
            Qcb = 0;
        end

        Qv = sum(h*(Tnow-Tc)*Avec.*vc); % Heat convection

        % now compute heat conduction from adjecent elements
        Tvec = zeros(1,6);  % get T of adject elements
        for i = 1:6
            if mynblist(i) > 0
                %    can this loop be changed to a vectorized command?

                Tvec(i) = T_now(mynblist(i));%Temp(tstep, mynblist(i));
            end
        end
        Tvec = Tvec - Tnow;
        kvec = know;
        Qc = sum(kvec.*Bvec.*cc.*Tvec);

        % laser
        if LaserSwitch && topelements(ia)
            Qi = total_Qi(ia); % irradiation
        else
            Qi = 0;
        end

        % new T
        %Temp(tstep+1,ia) = (Qc + Qcb - Qv)/C1/Cnow + Tnow;  %!! this C_p needs to be considered as T dependent.
        Tnew(ia) = (Qc + Qcb - Qv + Qi)/C1/Cnow + Tnow;  %new T of existing elements

    end % end loop elements


    %Temp(tstep+1,1:eNum) = Tnew;
    Temp_new(1:eNum) = Tnew;
    

    tstep = tstep + 1; % advance time
    eNum = eNum + 1; % deposite the next element
    cx = mypath(eNum,1);
    cy = mypath(eNum,2);
    cz = mypath(eNum,3);
    StateM(cx,cy,cz) = 1;
    %Temp(tstep,eNum) = Te;
    
    Temp_pre = Temp_new;
    Temp_pre(1,eNum) = Te;
    
    % need to write data!!!! here
    if mod(tstep-1,pr_interval) == 0
        fprintf(fsteps,'%d ',tstep);
        fprintf(fsteps,'%f ', Temp_pre);
        fprintf(fsteps, '\n');
        
        fwrite(fsteps2, tstep, intfmt);
        fwrite(fsteps2, Temp_pre, 'single');
    end
    % update nb_list due to the addition of current element
    nb_list(eNum,:) = 0; % init
    if cz ==1 % at first layer, conduction with bed
        nb_list(eNum,6) = -1;
    end
    % if not first element, start finding neighbours
    [nb_num, nb_face] = findnb(eNum, mypath);

    % now update nb list of current element and the nb list of nbs
    for i = 1:length(nb_face)
        nb_list(eNum, nb_face(i)) = nb_num(i);

        % 1 <--> 4, 2 <--> 5, 3 <--> 6
        nb_face_inv = nb_face(i) - 3 * abs(nb_face(i)-3.5)/(nb_face(i)-3.5);
        nb_list(nb_num(i), nb_face_inv) = eNum;
    end


end
tloope = toc;
fclose(fsteps); % close the data file
fclose(fsteps2);
% print a simulation summary to text log file
fname = getfname(); % auto-generate a file name
f1 = fopen(fname,'w');
if (f1 == -1) 
    fprintf("Error opening log file for writing data");
end
fprintf(f1,"*******************Simulation Done!*****************\n");
fprintf(f1,"Total run time\n");
fprintf(f1,"Space discretization dx = %fmm  dy = %fmm  dz = %fmm\n",dx*1e3,dy*1e3,dz*1e3);
fprintf(f1,"Time discretization dt = %fs\n",dt);
fprintf(f1,"Simulation box size Lx = %fmm  Ly = %fmm  Lz = %fmm\n",Lx*1e3,Ly*1e3,Lz*1e3);
fprintf(f1,"Simulation time: %d\n",dt*tstep);
fprintf(f1,"Total number of elements: %d\n",natoms);
tend = toc;
fprintf(f1,"*****run infomation*******\n");
%HWinfo = system('systeminfo');  % does not work for mac
fprintf(f1,"Arch: %s\n",computer);
fprintf(f1,"Total parallel workers: %d\n", numlabs);
fprintf(f1,"Main loop time = %f s\n", tloope-tloops);
fprintf(f1,"Total Run time = %f s\n", toc);

fclose(f1);

function fname=getfname()
% this function is used to get current time and generate a filename for
% output general simulation info.
ctime = clock;
fname = sprintf("Log_%4d%02d%02d%02d%02d.txt",ctime(1),ctime(2),ctime(3),ctime(4),ctime(5));
end

%{

fprintf("*******************Simulation Done!*****************")
fprintf("Total run time")
fprintf("Space discretization dx = %fmm  dy = %fmm  dz = %fmm\n",dx*1e3,dy*1e3,dz*1e3)
fprintf("Time discretization dt = %fs\n",dt)
fprintf("Simulation box size Lx = %fmm  Ly = %fmm  Lz = %fmm\n",Lx*1e3,Ly*1e3,Lz*1e3)
fprintf("Simulation time: %d\n",dt*tstep)
fprintf("Total number of elements: %d\n",natoms)
tend = toc;
fprintf("*****run infomation*******")
HWinfo = system('systeminfo');
fprintf("Total parallel workers: %d\n", numlabs)
fprintf("Main loop time = %f s\n", tloope-tloops)
fprintf("Total Run time = %f s\n", toc)

%}

