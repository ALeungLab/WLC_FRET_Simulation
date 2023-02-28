% WLC_FRET_simulation
% Adam Kenet, Mohsen Badiee
% modified from TJ Ha
% updated 2/4/2022


function WLC_FRET_simulation_AK(polymer,b0,const,computer,save)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Driver function that uses experimental data and simulates FRET to find persistence length for the polymer.
% To enter new data, create a new switch/case
%
% FUNCTIONS NEEDED:
%          1. optimize_AK
%             a. wlc_AK
%
% INPUTS:
%   polymer  -- [string] -- which switch-case to run (where is the experimental data) (eg. 'PAR', 'RNA')
%   b0       -- [int]    -- size of monomer in Angstroms (eg. 11.6)
%   const    -- [int]    -- step size for random walk in Angstroms (eg. 13.8)
%   computer -- [string] -- options: 'mac' or 'windows' (which computer are you using-need for saving the data)
%   save     -- [string] -- must be 'save' if you want to save the data; any other string will result in the data not being saved
%
% OUTPUTS (will save to files in the directory: root/Final_Data/polymer/b0/const)
%       eg: Desktop/Final_Data/PAR/11.6/13.8
%
%   DATA SAVED:   (NB: columns are decreasing salt concentration)
%       persistence_length -- [1xn array, n=number of salt conditions]
%                          -- persistence length for polymer
%       RMSD -- [1xn array, n=number of salt conditions]
%            -- root mean squared deviation for each salt concentration
%       E_sim -- [mxn matrix, m=number of lengths to simulate (from "range"), n=number of salt conditions]
%             -- simulated FRET values
%
%   saves persistence_length to file "persistence_length.csv"
%   saves RMSD to file "RMSD.csv"
%   saves E_sim to file "E_sim.csv"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

root = pwd; % get working directory

switch polymer
  
    case 'PAR'
        howmany=[7 12 17 22 27 35]; % sizes of polymer [-mer]

        R0=59;      % Forster distance [Angstrom]

        range=2:40; % -mer range to simulate

        numsalt=8; % number of salt conditions
        salt=1:numsalt; % salt condition

        pos_vals=10:10:180; % values of plArray to test (before fine tuning)

        % Experimental PAR data
        % ROW=Length
        % COLUMN=Salt
        E_exp = zeros(length(howmany),numsalt); % y-axis of plot
        E_exp(:,8)=[0.55  0.22  0.19  0.10  0.07  0.10]'; % 25 mM
        E_exp(:,7)=[0.56  0.24  0.19  0.13  0.08  0.11]'; % 50 mM
        E_exp(:,6)=[0.64  0.30  0.27  0.15  0.11  0.12]'; % 100 mM
        E_exp(:,5)=[0.69  0.44  0.28  0.23  0.15  0.17]'; % 200 mM
        E_exp(:,4)=[0.75  0.58  0.41  0.33  0.24  0.25]'; % 400 mM
        E_exp(:,3)=[0.80  0.68  0.55  0.44  0.39  0.34]'; % 800 mM
        E_exp(:,2)=[0.80  0.72  0.60  0.47  0.44  0.37]'; % 1000 mM
        E_exp(:,1)=[0.82  0.78  0.71  0.60  0.61  0.51]'; % 2000 mM
                  % 7mer  12mer 17mer 22mer 27mer 35mer

    case 'RNA'
        howmany=[20 30 40 50 70]; % sizes of polymer [-mer]

        R0=59;      % Forster distance [Angstrom]

        range=2:75; % -mer range to simulate

        numsalt=8; % number of salt conditions
        salt=1:numsalt; % salt condition
        
        pos_vals=10:10:180; % values of plArray to test (before fine tuning)
        
        % Experimental RNA data
        % ROW=Length
        % COLUMN=Salt
        E_exp = zeros(length(howmany),numsalt); % y-axis of plot
        E_exp(:,8)=[0.46  0.32  0.24  0.17  0.12]'; % 25 mM
        E_exp(:,7)=[0.51  0.30  0.24  0.19  0.12]'; % 50 mM
        E_exp(:,6)=[0.55  0.34  0.26  0.20  0.13]'; % 100 mM
        E_exp(:,5)=[0.57  0.39  0.30  0.23  0.14]'; % 200 mM
        E_exp(:,4)=[0.64  0.45  0.36  0.30  0.19]'; % 400 mM
        E_exp(:,3)=[0.69  0.52  0.43  0.36  0.25]'; % 800 mM
        E_exp(:,2)=[0.71  0.53  0.46  0.39  0.26]'; % 1000 mM
        E_exp(:,1)=[0.73  0.57  0.52  0.43  0.33]'; % 2000 mM
                  % 20mer 30mer 40mer 50mer 70mer

end


[plArray,RMSD,E_sim] = optimize_AK(E_exp,pos_vals,const,b0,R0,range,salt,howmany,numsalt);

% calculate persistence length from plArray
persistence_length = (13 + (plArray/4))/10; % [nm] % formula from TJ Ha

% save data
if strcmpi(save,'save') % only save the data if given the input
    % note: columns are decreasing salt concentration
    if strcmpi(computer,'mac')
        data_dir = strcat(root,'/Final_Data/',polymer,'/',num2str(b0),'/',num2str(const));
    elseif strcmpi(computer,'windows')
        data_dir = strcat(root,'\Final_Data\',polymer,'\',num2str(b0),'\',num2str(const));
    end
    mkdir(data_dir); cd(data_dir);
    writematrix(persistence_length,'persistence_length.csv')
    writematrix(RMSD,'RMSD.csv')
    writematrix(E_sim,'E_sim.csv')
    cd(root)
end
end

function [new_plArray, min_RMSD, E_sim_optimal] = optimize_AK(E_exp,pos_vals,const,b0,R0,range,salt,howmany,numsalt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that finds the optimal plArray (can convert to persistence length)
% given experimental FRET data and simulated FRET (calculated from WLC model).
%
% Simulates FRET data and finds best fit by minimizing RMSD from experimental FRET data
%
% To enter new data, create a new switch/case in the driver function "WLC_FRET_simulation_AK"
%
% NEEDED FOR FUNCTIONS:
%          1. WLC_FRET_simulation_AK (driver function)
%
% FUNCTIONS NEEDED:
%          1. wlc_AK
%
% INPUTS:
%   E_exp  -- [mxn matrix, m=number of lengths, n=number of salt conditions]
%          -- experimental FRET values
%          -- NB: column 1 should have highest salt concentration
%          -- NB: columns should be in order of decreasing salt concentration
%   pos_vals -- [1xp arrary] -- values of plArray to test before fine tuning
%            -- eg. 10:10:180
%   const -- [int] -- step size for random walk in Angstroms (eg. 13.8)
%   b0    -- [int] -- size of monomer in Angstroms (eg. 11.6)
%   R0    -- [int] -- % Forster distance in Angstroms (eg. 59)
%   range -- [1xq array] -- lengths of polymer to simulate [-mer]
%                        -- eg. 2:40
%   salt -- [1xn array, n=number of salt conditions to simulate] 
%        -- array of integers increasing by 1 for each salt condition (eg. 1:8)
%   howmany -- [1xr array] -- array containing the lengths of polymer for which there is experimental FRET data
%                          -- eg. [7 12 17 22 27 35]
%   numsalt  -- [int] -- number of salt conditions to simulate
%
% OUTPUTS 
%   new_plArray -- [1xn array, n=number of salt conditions]
%               -- optimal plArray that minimizes RMSD to experimental data
%               -- can convert to persistence length
%   min_RMSD -- [1xn array, n=number of salt conditions]
%            -- root mean squared deviation for each salt concentration given optimal plArray
%   E_sim_optimal -- [mxn array, m=number of lengths to simulate (from "range"), n=number of salt conditions]
%                 -- simulated FRET values given optimal plArray
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plArray = zeros(1,numsalt); % initialize empty array
difsq = zeros(1,length(howmany));  % initialize empty array
E_sim_all = zeros(max(range),numsalt); % simulated vals from WLC function
E_sim2save = zeros(max(range),numsalt); % temporarily hold values for tested plArray value
E_sim_store = zeros(max(range),max(pos_vals)+1); % simulated values for optimal plArray value

counter2=1;
for i=salt % loop for each column in plArray
    counter=1; % initial tuning counter
    counter_FT=1; % fine tuning counter
    for j=pos_vals % loop to change value in one column of plArray
        disp(j)
        plArray_test = j; % changes ith value of plArray
        E_sim_all = wlc_AK(const,plArray_test,b0,R0,range,i,numsalt); % calculates Esim 

        E_sim=[]; % initiallize array
        for x=1:length(howmany) % for each experimental length
            val = E_sim_all(howmany(x),i);
            E_sim = [E_sim, val]; % append to array
            difsq(x) = (E_exp(x,i)-E_sim(x)).^2; % difference squared
        end
        rmsd = sqrt(mean(difsq)); % RMSD values
        attempted_vals(counter)=j; % attempted plArray values % matrix of j's % only contains j values
        rmsd_vals(counter)=rmsd; % matrix of rmsd % only contains rmsd values
        counter=counter+1;
    end
    [~,I_first]=min(rmsd_vals); % want index of min rmsd

    % fine tuning (+- 10 from initial min)

    for k=attempted_vals(I_first)-10:attempted_vals(I_first)+10 % find min from +- 10 of initial min
        disp(k)
        plArray_test = k; % changes ith value of plArray
        E_sim_all = wlc_AK(const,plArray_test,b0,R0,range,i,numsalt); % calculates Esim

        E_sim=[]; % initiallize array % simulated values at points where have experimental data
        for x=1:length(howmany) % for each experimental length
            val = E_sim_all(howmany(x),i);
            E_sim = [E_sim, val]; % append to array
            difsq(x) = (E_exp(x,i)-E_sim(x)).^2; % difference squared
        end
        E_sim_store(:,k+1)=E_sim_all(:,i); % temporarily store E_sim_all value % plus one because first val is zero and indexing starts at 1
        rmsd_FT = sqrt(mean(difsq));
        attempted_vals_FT(counter_FT)=k; % attempted plArray values % matrix of k's % only contains k values
        rmsd_vals_FT(counter_FT)=rmsd_FT; % matrix of rmsd % only contains rmsd values
        counter_FT=counter_FT+1;
        disp(attempted_vals_FT)
        disp(rmsd_vals_FT)
    end
    [minRMSD,I]=min(rmsd_vals_FT); % want index of min rmsd

    RMSDarray(counter2) = minRMSD;
    counter2=counter2+1;
    best_plArray_val = attempted_vals_FT(I); % matrix of j's (at index of min rmsd) % tells us the best value for plArray(i)
    plArray(i) = best_plArray_val;
    disp(plArray)
    disp(RMSDarray)

    E_sim2save(:,i)=E_sim_store(:,best_plArray_val+1); % keep the simulated values for the best plArray value

end % end salt loop
new_plArray = plArray;
min_RMSD = RMSDarray;
E_sim_optimal = E_sim2save;
end % end function

function [E_sim] = wlc_AK(const,plArray,b0,R0,range,salt,numsalt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that simulates the Worm Like Chain model and generates FRET values
%
% (in "wlc_AK"): is ran for each value of plArray that is being tested 
%
% modified from TJ Ha
%
% NEEDED FOR FUNCTIONS:
%          1. WLC_FRET_simulation_AK (driver function)
%          2. wlc_AK (optimization function)
%
% FUNCTIONS NEEDED:
%          none
%
% INPUTS:
%   const -- [int] -- step size for random walk in Angstroms (eg. 13.8)
%   plArray -- [] -- 
%   b0    -- [int] -- size of monomer in Angstroms (eg. 11.6)
%   R0    -- [int] -- % Forster distance in Angstroms (eg. 59)
%   range -- [1xq array] -- lengths of polymer to simulate [-mer]
%                        -- eg. 2:40
%   salt -- [1xn array, n=number of salt conditions to simulate] 
%        -- array of integers increasing by 1 for each salt condition (eg. 1:8)
%   numsalt  -- [int] -- number of salt conditions to simulate (eg. 8)
%
% OUTPUTS 
%   E_sim -- [mxn array, m=number of lengths to simulate (from "range"), n=number of salt conditions]
%         -- simulated FRET values given plArray
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nP=numsalt; %  nP = 8 salt conditions
nL=max(range)-min(range)+1;
finalArray=zeros(nL,nP);
inc=const*10^-10;
for m=salt % m = current salt condition

    for k=range % for each length

        L = ((k*b0)+12.2)*10^-10; % corrected formula for cy3(6.3A) and AMP(5.9A) on ends

        % calculate persistence length
        % note: this formula is "arbitrary" and is to ensure persistence length is positive 
        if length(plArray)==1 % for one salt condition
            lp = (13+0.25*plArray)*10^-10; % formula from TJ Ha
        else % for multiple salt conditions
            lp = (13+0.25*plArray(m))*10^-10; % formula from TJ Ha
        end

        t = L/lp;
        R0_conv=R0*10^-10;
        
        n = 4*(3 *t /4)^(3/2)*exp(3 *t /4)/(pi^(3/2)*(4 + 12/(3* t /4) + 15/(3* t /4)^2));
        N=1000;
        
        rArray=(1:N)/(N+1);
        
        probArray=zeros(N,1);
        sumArray=zeros(N,1);
        inverseArray=zeros(N,1);
        for i=1:N
            r=rArray(i);
            probArray(i)=4*pi*n*r^2/(1 - r^2)^(9/2)*exp(-3*t/4/(1 - r^2));
        end
        
        sumArray(1)=probArray(1);
        for i=2:N
            sumArray(i)=sumArray(i-1)+probArray(i);
        end
        sumArray=sumArray/N;
        
        j=1;
        for i=1:N
            while sumArray(j) < i/N
                j=j+1;
            end
            inverseArray(i)=rArray(j);
        end
        
        bigN=10000;
        
        transfer=0;
        for i=1:bigN
            alive=1;
            oldr=inverseArray(ceil(rand*N))*L;
            while alive > 0
                if rand <0.5
                    newr=oldr+inc;
                else
                    newr=oldr-inc;
                end
                
                if newr < L && newr > 0
                    
                    if rand < probArray(ceil(newr/L*N))/probArray(ceil(oldr/L*N))
                        oldr=newr;
                    end
                end
                
                rN=rand;
                if rN < 0.01
                    alive=0;
                else
                    if rN < 0.01+0.01/((oldr/R0_conv)^6)
                        alive=0;
                        transfer=transfer+1;
                    end
                end
            end
        end %closing bigN loop
        finalArray(k,m)=transfer/bigN;
    end %closing k loop
     disp(m);
end %closing m loop
E_sim=finalArray;
disp('SIMULATION DONE')
end