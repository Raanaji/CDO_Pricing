function [PV_def, PV_premium] = cash_flow(expiry,def_time,rec,zc_rate,capital,C,D)

% This function computes the value of the premoum (DV01, spread still to be
% determined) and the default legs of a CDO 
% Inputs are:
% expiry: CDO maturity
% def_time: simulated default time
% rec: recovery
% zc_rate: constant ZC rate
% capital: nominal amount of each reference entity
% C : attachment point
% D : detachment point

PV_def = 0;
PV_premium = 0;
num = size(def_time,2);
loss = zeros(num,1); % LGD for each credit
tot_loss = 0; % cumulative portfolio loss
periodic_loss = zeros(expiry,1); % stores the accumulated loss at each
                                 % payment date
out_capital = zeros(expiry,1);   % outstanding tranche capital at each
                                 % payment date
fee = zeros(expiry,1);
total_fee = 0;
indicator = 0;
c = 0;

% calculate the total loss in the k-th simulation
for i = 1:num
    if def_time(1,i)<expiry % if the simulated default time for the generic 
                            % credit is < than the CDO maturity, there is 
                            % a loss
       loss(i) = (1-rec)*capital;
       tot_loss = tot_loss + loss(i); % sum the individual losses
    end
end

%% DEFAULT LEG SIMULATION %%

% if the loss is below the lower treshold C there's no default payment

if tot_loss<C
    PV_def = 0;
    
% if the loss is above C and below D there's a default payment
elseif tot_loss>C & tot_loss<D
        for i=1:num
            if def_time(1,i)<expiry
                indicator = indicator + loss(i); % we memorize the cumulative
                                                 % losses
                if c==0
                    disc_fact_def = 0;
                    r = zc_rate;
                    disc_fact_def = (1+r)^(-def_time(1,i)); % discount factor
                                                             % at default
                    PV_def = PV_def + (indicator-C)*disc_fact_def; % only the 
                      % loss exceeding C is absorbed (look at the footnote
                      % 32 in the simulation algorithm p. 51 of my paper)
                    c = 1;
                else
                    dis_fact_def = 0;
                    r = zc_rate;
                    disc_fact_def = (1+r)^(-def_time(1,i));
                    PV_def = PV_def + loss(i)*disc_fact_def;
                end 
            end
        end
    
    
    % if portfolio losses are exceeding D, the C-D% tranche only absorb
    % losses upto D%
    elseif tot_loss>D
    for i = 1:num
        if def_time(1,i)<expiry
            indicator = indicator+loss(i);
            if indicator>C & indicator<D % check if losses are in the C-D% range
               
                if c==0
                    disc_fact_def = 0;
                    r = zc_rate;
                    disc_fact_def = (1+r)^(-def_time(1,i));
                    PV_def = PV_def + (indicator-C)*disc_fact_def;
                    c=1;
                else
                    disc_fact_def = 0;
                    r = zc_rate;
                    disc_fact_def = (1+r)^(-def_time(1,i));
                    PV_def = PV_def + loss(i)*disc_fact_def;
                    
                end
            elseif indicator>D % look at footnote 33 p. 52 of my paper
                if c==1
                    disc_fact_def = 0;
                    r = zc_rate;
                    disc_fact_def = (1+r)^(-def_time(1,i));
                    absorbed_loss = D - (indicator-loss(i))*disc_fact_def;
                    % look at footnote 33 p . 52 of my paper
                    PV_def = PV_def+(absorbed_loss*disc_fact_def);
                    c=2;
                end
            end
        end
    end
end 


%% PREMIUM LEG SIMULATION %%

for i=1:expiry
    periodic_loss(i)=0;
    for j=1:num
        if def_time(1,j)<1
            % calculated the accumulated portfolio losses
            periodic_loss(i) = periodic_loss(i)+(1-rec)*capital;
        end
    end 
    out_capital(i) = min(max(D-periodic_loss(i),0),D-C); % outstanding capital at each payment date
    
    fee(i) = ((1+zc_rate)^(-i))*out_capital(i);
    PV_premium = PV_premium + fee(i); % DV01
    
end 


        
                    
                    
                  
            

