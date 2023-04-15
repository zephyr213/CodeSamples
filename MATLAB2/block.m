classdef block < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    % 10/12 replaced obj name as bk

    properties
        ntot = 0; % expected total elements.
        eNum = 0; % count of current deposited elements
        Temp_pre = nan; % temperature of all elements from (N-1) step
        Temp_new = nan; % temperature of all elements from N step
        Cxyz;  % coordinates of all elements
        nblist = nan; % neighbour list. 6 x ntot.  1- 2- 3- 4- 5- 6-

        newelist = [];   % list of newly added elements, pending nb list update, row vector

        % temp_record related
        Temp_record = nan; % variable used to temperarily store the Temp of each elements
        maxRecordStorage = 2; % GB, max size of Temp_record
        maxRecordRows; % max number of steps saved in Temp_record
        recordRowCounter = 1; % count the row to be written to in Temp_record
        recordInterval = 1; % step intervals to record Temperature.

        % 
        istep = 0; % current simulation step
        nsteps = 100; % total number of simulation steps

        % output file related
        fname; % file name
        fid;    % file id

        % switch for physics:
        switchPhysics = struct('Laser', false, ... % for laser irradiation
            'Qc', true, ... % for heat conduction with neighbours
            'Qv', true, ... % for heat convection
            'Qcb', true);   % for heat conduction with bed

        % switch for other features
        switchFeatures = struct(...
            'checkDuplicate', false ... % check if newly deposited elements overlaps with existing ones
            );

        % element dimensions 1x3 row vector
        % coordinates of last deposited element?


    end

    properties (SetAccess=private)
        % simulation parameters
        % need to write a function to se these simulation parameters
        dxyz = struct('dx',[],'dy',[],'dz',[]);
        dt = []; % s, timestep size


        % processing parameters
        ppara = struct('Tb',[],'Tc',[],'Te',[],'h',[],'v',[]);

        % materials parameters
        mpara = struct('Rho',[],'kx',[],'ky',[],'kz',[]);

        % laser parameters
        lpara = struct('S',[],'R',[],'P',[],'G',[],'ab',[]);

        % constants used for computation
        % needs to call to update the constants after all parameters are
        % entered.
        const_C = struct('C1',[],'C2',[]);
        const_vec = struct('Avec0',[], 'dvec',[]);
    end

    methods
        function bk = block(N)
            %block Construct an instance of this class
            %   Initialize important variables
            bk.ntot = N;
            bk.Temp_new = nan(1,N);
            bk.Temp_pre = bk.Temp_new;
            bk.Cxyz = nan(3,N);
            bk.nblist = nan(6,N);

            % determine maxstorage value
            bk.maxRecordRows = floor(bk.maxRecordStorage * 1e9/(N+1)/8);
            bk.Temp_record = nan(bk.maxRecordRows,N+1); % first column for steps

            % get filename and open file
            bk.fname = bk.getfname;
            bk.fid = fopen(bk.fname,"w");
            if bk.fid == -1
                fprintf("Error opening log file for writing data!!")
            end

            % need to call other function to update the Mpara, Ppara, Spara,Lpara
            % currently just use default values as a test
 

        end

       % function bk = updateMpara(varargin)
       function bk = updateMpara(bk,varargin)

            % update materials parameters
            % currently set to assign constant to the parameters
            if nargin == 1
                bk.mpara.Rho = 1125; % kg/m3
                bk.mpara.kx = 0.17; bk.mpara.ky = 0.17; bk.mpara.kz = 0.17; % W/mK, thermal conductivity along x,y,z
            end

       end
        function bk = updatePpara(bk,varargin)
            % update processing parameters
            % currently set to assign constant to the parameters
            % processing parameters
            if nargin == 1
                c2k = 273.15;
                bk.ppara.h = 100; % W/m2K, coefficient of convection
                bk.ppara.Tb = 90 + c2k; % bed T
                bk.ppara.Tc = 40 + c2k; % chamber T
                bk.ppara.Te = 240 + c2k; % extrusion T
                bk.ppara.v = 50e-3; % extrusion speed 50 mm/s taken from Garon's notes.
            end
        end
        function bk = updateSpara(bk,varargin)
            % update simulation parameters
            % varargin: dx,dy,dz,dt,nsteps
            if nargin == 1 || nargin ~=6
                % some default values
                bk.dxyz.dx = 0.01;
                bk.dxyz.dy = 0.01;
                bk.dxyz.dz = 0.01;
                bk.dt = 0.01;
                bk.nsteps = 4;
            else
                bk.dxyz.dx = varargin{2};
                bk.dxyz.dy = varargin{3};
                bk.dxyz.dz = varargin{4};
                bk.dt = varargin{5};
                bk.nsteps = varargin{6};
            end
        end

        function bk = updateLpara(bk,lswitch,varargin)
            % update laser parameters
            % some default values
            fprintf("%d\n",nargin)
            bk.switchPhysics.Laser = lswitch;
            if lswitch && nargin == 2 
                v = bk.ppara.v; 
                bk.lpara.S = 11*v;  % laser speed
                bk.lpara.R = 4e-4/2;    % m, laser focal radius, data from experiment d=403 micrometers
                bk.lpara.P = 1.5; % W, laser power
                bk.lpara.G = bk.lpara.P/pi/bk.lpara.R^2;    % laser irradiation W/m2
                bk.lpara.ab = 0.2  ;   % laser abosorption
                % set 0.2 for 10.60 micrometer IR laser as a test
            end
        end

        function bk = computeConstants(bk)
            % copy properties to local function
            dx = bk.dxyz.dx;
            dy = bk.dxyz.dy;
            dz = bk.dxyz.dz;
            dt = bk.dt; %#ok<PROP> 
            Rho = bk.mpara.Rho;
            kz = bk.mpara.kz;
            % face vector [Ax,Ay,Az,Ax,Ay,Az]
            Avec0 = [dy*dz, dx*dz, dx*dy];
            Avec0 = [Avec0, Avec0]';
            %kvec = [kx,ky,kz,kx,ky,kz];
            dvec = [dx,dy,dz,dx,dy,dz]';
            %Bvec = Avec./dvec;
            C1 = Rho * dx * dy * dz/dt; %#ok<PROP> 
            C2 = kz*dx*dy/dz;

            bk.const_C.C1 = C1;
            bk.const_C.C2 = C2;
            bk.const_vec.Avec0 = Avec0;
            bk.const_vec.dvec = dvec;
        end

        function bk = addElements(bk,n,x,y,z,T)
            %add element(s) to a block
            %   n - number of elements to be added
            %   x,y,z - coordinates of elements (needs to be row vectors if
            %   n>1)
            %   T - temperature of new elements (needs to be row vectors if
            %   n>1)

            % optional check duplicate elements
            if bk.switchFeatures.checkDuplicate
                if existDuplicate(bk.Cxyz,[x;y;z])
                    fprintf("Warning: new deposition at step: %d\n overlap with existing ones!!\n", bk.istep)
                    fprintf("WARNING: No new elements deposited!")
                    return % make sure this will exit the function
                end
            end

            eNum_old = bk.eNum;
            if any([numel(x), numel(y), numel(z), numel(T)] ~= n)
                fprintf("WARNING: x,y,z sizes don't match with element number entered")
                fprintf("WARNING: No new elements deposited!")
                return;
            end

            % Now start adding new elements:
            bk.eNum = bk.eNum + n;    % increase deposited elements
            bk.Temp_pre(eNum_old+1:bk.eNum) = T;  % add the T from new elements
            bk.Cxyz(:,eNum_old+1:bk.eNum) = [x;y;z];  % add the T from new elements

            bk.newelist = [bk.newelist, eNum_old+1:bk.eNum]; % update new elements list (will set to empty when updateNBlist is called.)
            if bk.eNum > bk.ntot   % case when newly deposited elements increase total above preset ntot.
                ndiff = bk.eNum -bk.ntot;
                fprintf("WARNING: deposited elements greater than preset total, diff=%d\n",ndiff);
                fprintf("WARNING: block size automatically increased from %d to %d\n",bk.ntot, bk.ntot+n);
                bk.Temp_pre = [bk.Temp_pre, nan(1,ndiff)];
                bk.ntot = bk.ntot + ndiff;
            end
        end

        function isduplicate = existDuplicate(oldXYZ, newXYZ)
            % check if any elements in new set overlap with existing ones
            % oldXYZ -- 3 x n matrix, coordinates of existing elements
            % newXYZ -- 3 x m matrix, coordinates of new elements
            isduplicate = false;
            for e = newXYZ
                if any(all(oldXYZ == e))
                    isduplicate = true;
                    return;
                end
            end

        end

        function bk = updateNBlist(bk,varargin)
            %update nblist
            %   varargin - by default, only update nblist of newly added
            %   elements and their neighbors
            if nargin == 2 && varargin{2} == "all"
                e2update = 1:bk.eNum;
            else
                e2update = bk.newelist;
            end

            % start update nblist
            bk.nblist(:,e2update) = 0; % initialize as 0 (0=convection to air by default)
            %bk.nblist(6,e2update(bk.Cxyz(3,e2update) == 1)) = -1; % bottom layer, conduction with bed

            % define a matrix used in for loop
            mat1 = [1 0 0; 0 1 0; 0 0 1];
            mat1 = [mat1; -mat1];

            % possible to get rid of FOR loop using array operations?
            for i = e2update  % go over each "new" element
                % bottom layer, conduction with bed. - check with above
                % vectorized operation, will it be more efficient?
                if bk.Cxyz(3,i) == 1
                    bk.nblist(6,i) = -1;
                end

                dis = bk.Cxyz - bk.Cxyz(:, i);
                dis_abs = abs(dis);
                % only consider immediate neighbors here. If need to
                % consider additinal neighbors, change the way nl is
                % determined
                nl = (sum(dis_abs,1) == 1); % logical vector to select neighbours
                nl_dis = dis(:, nl);    % distance matrix of neighbors

                tmpn = 1:bk.eNum;
                nb_num = tmpn(nl); % List of neighbours' indices

                % now find the face of the neighbors

                c = logical(mat1*nl_dis == 1);
                tmp_contact = repmat((1:6)',1,size(c,2));
                nb_face = tmp_contact(c)'; % List of contact face number, rv

                % now update nb list of current element
                bk.nblist(nb_face,i) = nb_num';

                % and the nb list of nbs
                nb_face_inv = nb_face - 3 * abs(nb_face-3.5)./(nb_face-3.5);
                ind = sub2ind(size(bk.nblist),nb_face_inv,nb_num);
                bk.nblist(ind) = i;
            end
            bk.newelist = []; % reset updatelist

        end

        function bk = evolveT2(bk)
            % design this function to evolve T of certain elements
            % for future development, can just modify the elist variable
            % (e.g. set cutoff distance to allow elements far away from the newly deposited
            % to be treated as constant T or update infrequently for better
            % efficiency

            % check if nblist updated
            if ~isempty(bk.newelist)
                fprintf("Warning: neighbor list not updated after new elements deposition at Step: %d\n", bk.istep);
            end

            % record initial conditions?
            if bk.istep == 0
                bk.Temp_record(bk.recordRowCounter,:) = [bk.istep, bk.Temp_pre];
                bk.recordRowCounter = bk.recordRowCounter + 1;
            end


            %evolve system T by one step
            bk.istep = bk.istep + 1;

            if bk.istep > bk.nsteps % if evovle T beyond nsteps
                fprintf("!!!Going beyond preset total steps!Automatically increase nsteps by 1!!!\n")
                bk.nsteps = bk.nsteps + 1;
                % automatically increase steps? Or stop?
            end


            % determine the elements to be updated
            %   first select all the eNum elements
            elist = ~isnan(bk.Temp_pre);
            ne = sum(elist);
            % should use ne and elist below.
            % ne as # of elements to be updated. elist as the indices
            % of elements need to be treated.

            Tpre = bk.Temp_pre(elist);
            Tpre6 = repmat(Tpre,6,1);
            mynblist = bk.nblist(:,elist); % nb list of current elements
            % Tnew = zeros(1,ne); % initialize

            % get k and C at current Temperature from interpolation
            % needs double check this part about the fitting
            [kvec,Cnow] = bk.fitkc(Tpre);

            % P1: calculate Qcb - conduction from bed ************
            Qcb = zeros(1,ne);
            Tb = bk.ppara.Tb; % ia
            C2 = bk.const_C.C2;

            % calculation
            bottomE = mynblist(6,:) == -1; % select bottom elements
            if bk.switchPhysics.Qcb
                Qcb(bottomE) = C2*(Tb - Tpre(bottomE));
            else
                Qcb(bottomE) = 0;
            end


            % P2: calculate Qc - conduction from neighbors ************
            % temp variables? or initialize variables
            Tvec = zeros(6,ne);

            % temp variables, needs to set as properties later.
            %dx = 1; dy = 2; dz = 3; % temp values as test

            %Avec = [dy*dz, dx*dz, dx*dy];
            %Avec = [Avec, Avec]';
            Avec0 = bk.const_vec.Avec0;
            Avec = repmat(Avec0,1,ne);  % match the dimensions of (6,ne)
            %dvec = [dx,dy,dz,dx,dy,dz]';
            if bk.switchPhysics.Qc
                dvec = repmat(bk.const_vec.dvec,1,ne);
                Bvec = (Avec./dvec);

                % mynblist > 0, select all faces of all elements with adjacent
                % elements
                cfaces = mynblist > 0;
                Tvec(cfaces) = Tpre(mynblist(cfaces));
                Tvec = (Tvec - Tpre6).*cfaces;

                % now Qc
                Qc = sum(kvec .* Bvec .* Tvec);
                % note here kvec has 1 row, while Bvec and Tvec have 6 rows.
                % MATLAB will auto expand kvec to match the size
            else
                Qc = zeros(1,ne);
            end

            % P3: calculate Qv - convection from env ************
            % variables
            if bk.switchPhysics.Qv
                Tc = bk.ppara.Tc; % env Temp
                h = bk.ppara.h; % convection parameter

                vfaces = mynblist == 0;
                Qv = sum(h .* (Tpre6 - Tc) .* Avec .* vfaces);
            else
                Qv = zeros(1,ne);
            end


            % P4: calculate Qi - laser irradiation ************
            % assume laser speed is faster than deposition head speed?
            if bk.switchPhysics.Laser
                % calulate the irradiation from laser assisted heating
                Qi = zeros(1,ne);

                % identify all elements affected by laser.
                % currently assume only the top layer elements
                
                topE = mynblist(3,:) == 0; % select top elements

                % now determine the paths of the laser and those top
                % elements being affected

                % needs properties to record laser coordinates, 
                % needs a function to govern the movement of laser head
                % then determine how many times each top element is
                % irradiated
                % try to use vectorization for better efficiency

            else
                Qi = zeros(1,ne);  % no laser heating
            end


            % Pfinal: add together.
            C1 = bk.const_C.C1;
            % evolve the temperature
            Tnew = (Qc + Qcb - Qv + Qi)./C1./Cnow + Tpre;
            % update property
            bk.Temp_new(elist) = Tnew;

            % Record Temperature 
            bk.Temp_record(bk.recordRowCounter,:) = [bk.istep, bk.Temp_new];
                % if additional values need saving, can be added to the
                % line above. Need to make sure the size of Temp_record is
                % correct
            bk.recordRowCounter = bk.recordRowCounter + 1;

            if bk.recordRowCounter > bk.maxRecordRows
                % call a function to output Temp_records to file if
                % Temp_record is at max size
                bk.writedata(bk.fid, bk.Temp_record);
                bk.Temp_record = nan(size(bk.Temp_record));
                bk.recordRowCounter = 1;
            end

            % update Temp_pre 
            bk.Temp_pre = bk.Temp_new;  
        end


        function endsimulation(bk)
            % function to close files if simulation ends.
            % also dump remaining data in Temp_record to file
            bk.writedata(bk.fid, bk.Temp_record);
            closestate = fclose(bk.fid);
            if closestate == -1
                fprintf("Error closing output file!!");
            end
        end
    end

    methods(Static)
        function [knew,Cnew] = fitkc(T)
            % find temperature dependent k and Cp
            % need to find a more efficient way to do the fit. Table
            % lookup?


            % ABS data https://www.sciencedirect.com/science/article/pii/S0264127519308470?via%3Dihub#f0010
            T0 = 0:50:250;
            T0 = T0 + 273.15;
            k0 = [.23,.25,.28,.29,.31,.33];
            C0 = [780,1040,1490,1710,1865,2020];

            % use linear interpolation as a test for now
            knew = interp1(T0,k0,T);
            Cnew = interp1(T0,C0,T);
        end

        function writedata(fid, Tmat)
            % function to write Temp data to a text file
            % h_file - file handle
            % Tmat - temperature matrix to be written to a file
            datamat = Tmat(~isnan(Tmat(:,2)),:);  % choose existing data
            fprintf(fid, ['%d\t' repmat('%10.5f\t', 1, size(datamat,2)-1) '\n'], datamat');
        end

        function fname=getfname()
            % this function is used to get current time and generate a filename for
            % output general simulation info.
            ctime = clock;
            fname = sprintf("Log_%4d%02d%02d%02d%02d.txt",ctime(1),ctime(2),ctime(3),ctime(4),ctime(5));
        end
    end
end