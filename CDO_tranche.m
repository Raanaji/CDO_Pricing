function [price_eq,price_mezz,price_senior] = CDO_tranche(ref_ent,T,k)

% this function computes the fair price of a synthetic CDO tranched in
% three portions, 0-3%, 3-14%, and 14-100%. Spread for each reference
% entity is fixed @ 150bps. We let recoveries and pairwise correlation
% changing in order to produce a chart of CDO breakeven spread vs.
% correlation recovery

% Inputs:

% ref_ent: # of reference netities in the pool
% T:       # maturity of the deal
% k:       # of simulation


tic
% Initialize a vector of zeros for each of the five tranches

CDO_P1 = zeros(5,5);
CDO_P2 = zeros(5,5);
CDO_P3 = zeros(5,5);
CDO_P4 = zeros(5,5);
CDO_P5 = zeros(5,5);

% Initialize a vector of zeros for different recoveries

Recovery = zeros(ref_ent,5);
hazard = zeros(5,1);

% Set the obligors' spread to 150bps p.a.

spread = 150/10000;

% Since the CDS term structure is flat, the relationship hazard rate =
% spread/(1-recovery) holds. Threrfore, fore each recovery we calculate the
% corresponding hazard rate.

for rec_cycle = 1:5
    Recovery(:,rec_cycle) = (.2*rec_cycle) - .2; % recovery changes from 0% to 80%
    hazard(rec_cycle) = spread/(1-Recovery(1,rec_cycle));
end

ZC = 0.05; % constant interest rate

Amount = zeros(ref_ent,1); % vector of notional amount for each credit
Amount (:)= 100; % each credit has a notional amount of 100 units
C = zeros(5,1); % we fix five attachment points: 0%, 3%, 6%, 12%, 22%
D = zeros(5,1); % we fix detachment points: 3%, 6%, 12%, 22%, 100%
C(1) = (0/100)*sum(Amount);
D(1) = (3/100)*sum(Amount);
C(2) = (3/100)*sum(Amount);
D(2) = (6/100)*sum(Amount);
C(3) = (6/100)*sum(Amount);
D(3) = (12/100)*sum(Amount);
C(4) = (12/100)*sum(Amount);
D(4) = (22/100)*sum(Amount);
C(5) = (22/100)*sum(Amount);
D(5) = (100/100)*sum(Amount);

time = zeros(ref_ent,1);
index = zeros(ref_ent,1);

R = [0:.2:0.8]; % constant pairwise correlation ranges from 0% to 80%
% start the correlation loop

for R_cycle = 1:5
    for xx = 1:ref_ent
        for yy = 1:xx
            if xx==yy
                corr(xx,yy) = 1;
            else
                corr(xx,yy) = R(R_cycle); % populate the correlation matrix
                corr(yy,xx) = R(R_cycle);
            end
        end
    end 
    
    def_t = gaussian_time(corr,k,ref_ent); % generate pseudodefault times
                                           % with gaussian copula and
                                           % constant hazard rate 
                                           
S_fees = zeros(5,5); % dummy variable for memorizing the simulated payment leg 
S_default = zeros(5,5); % dummy variable for memorizing the simulated default leg
M_fees = zeros(5,5); % variable which meorize the payment leg for each loop of recovery & correlation
M_default = zeros(5,5); % variable which memorize the payment leg for each loop of recovery & correlation

for n = 1:k % start the simulation loop
    for rec_cycle = 1:5 % start the recovery loop
        [time,index] = sort(def_t(n,:)); % sort the pseudo vector of default times
        tau = [time./hazard(rec_cycle);index]; % generate the vector of default time by dividing the 
                                               % corresponding hazard rate
                                               
        for u = 1:5 % start the loop for each of the three tranches 
        recovery = 0;
        fees = 0;
        % calculate the premium and default leg
        [default,fees] = cash_flow(T,tau,Recovery(1,rec_cycle),ZC,Amount(1),C(u),D(u));
        S_fees(rec_cycle,u) = S_fees(rec_cycle,u) + fees;
        S_default(rec_cycle,u) = S_default(rec_cycle,u)+default;
        end
    end 
end 

for u = 1:5
    for rec_cycle = 1:5;
    M_fees(rec_cycle,u)=S_fees(rec_cycle,u)/k; % average DV01
    M_default(rec_cycle,u) = S_default(rec_cycle,u)/k; % average default leg
    end
end

for rec_cycle = 1:5
    CDO_P1(R_cycle,rec_cycle) = (M_default(rec_cycle,1)/(M_fees(rec_cycle,1)))*10000; % B\E spread 
                                                                            % for the 0-3% tranche
    CDO_P2(R_cycle,rec_cycle) = (M_default(rec_cycle,2)/(M_fees(rec_cycle,2)))*10000; % B\E spread 
                                                                            % for the 3-6% tranche
    CDO_P3(R_cycle,rec_cycle) = (M_default(rec_cycle,3)/(M_fees(rec_cycle,3)))*10000; % B\E spread 
                                                                           % for the 6-12% tranche
    CDO_P4(R_cycle,rec_cycle) = (M_default(rec_cycle,4)/(M_fees(rec_cycle,4)))*10000; % B\E spread 
                                                                           % for the 12-22% tranche
    CDO_P5(R_cycle,rec_cycle) = (M_default(rec_cycle,5)/(M_fees(rec_cycle,5)))*10000; % B\E spread 
                                                                         % for the 22-100% tranche                                                                        % for the 14-100% tranche                                                                            
end
end

price_eq = CDO_P1;
price_mezz1 = CDO_P2;
price_mezz2 = CDO_P3;
price_senior = CDO_P4;
price_SuperSenior = CDO_P5;

figure(1)
surf((0:.2:.8),R,price_eq);
title('Equity Tranche (0%-3%)')
xlabel('Recovery');
ylabel('Correlation');
zlabel('Tranche spread (bps per annum)');

figure(2)
surf((0:.2:.8),R,price_mezz1);
title('Mezzanine Tranche 1 (3%-6%)')
xlabel('Recovery');
ylabel('Correlation');
zlabel('Tranche spread (bps per annum)');

figure(3)
surf((0:.2:.8),R,price_mezz2);
title('Mezzanine Tranche 2 (6%-12%)')
xlabel('Recovery');
ylabel('Correlation');
zlabel('Tranche spread (bps per annum)');

figure(4)
surf((0:.2:.8),R,price_senior);
title('Senior Tranche (12%-22%)')
xlabel('Recovery');
ylabel('Correlation');
zlabel('Tranche spread (bps per annum)');

figure(5)
surf((0:.2:.8),R,price_SuperSenior);
title('Super Senior Tranche (22%-100%)')
xlabel('Recovery');
ylabel('Correlation');
zlabel('Tranche spread (bps per annum)');

toc
                                                                            
                                                                            
   
    

                                           
                                         
                                           
                                        
                                    
                                          
