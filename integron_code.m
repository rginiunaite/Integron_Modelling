clear all
%% Integron evolutionary dynamics

%% parameters of the model 

n = 3; % number of different gene cassettes
k = 3; % max number of cassetes in integron
K = 10^9; % carrying capacity
n_0 = 0.1; % baseline death rate
n_i = 0.001; % increase in death rate caused by functional integrase
n_s = 0.3; % increase in death rate caused by stressor in absence of resistance gene
ro = 0.001; % rate of integrase-mediated gene reshuffling
theta = 0.5; % probability of re-inserting an excised gene
beta = 0.5; % how fast gene expression declines with increasing distance from promoter
gamma = 3; % parameter that determines how expression level of a resistance gene affects death rate
sigma_m = 0.2; % average fraction of time that a stressor is present
sigma_v = 0.01; % average rate of switches between presence and absence
mu = 10^-5; % mutation rate from functional ro non-functional
phi = 0.3; % threshold level when integrase becomes active
tau = 10^-16; % rate of HGT
S = [0;0;0]; % initial stress
time_step = 0;
final = [];
t_final = [];


count=1;


%%
% number of integron genotypes
N = n*(n-1)*(n-2)+n*(n-1)+n + 1;

% initialise bacteria with different genotypes and functional integrase

X = zeros(N,1);

% matrix that shows which cassetes are present in all three different
% positions


X_cas = [0,0,0;1,0,0;2,0,0;3,0,0;1,2,0;1,3,0;...
    2,1,0;2,3,0;3,1,0;3,2,0;1,2,3;1,3,2;2,1,3;2,3,1;3,1,2;3,2,1;0,0,0;1,0,0;2,0,0;3,0,0;1,2,0;1,3,0;...
    2,1,0;2,3,0;3,1,0;3,2,0;1,2,3;1,3,2;2,1,3;2,3,1;3,1,2;3,2,1];





% stress transition matrix
 
M = sigma_v / 2 * [1 - (1/(1-sigma_m)),1/(1-sigma_m);1/sigma_m, 1 - 1/sigma_m];


    
    
% expression level at different positions

for i =1:k
    E(i) = exp(-beta*(i-1));
end

% total gene expression of cassette j

% E_total=zeros(2*N,n);
% 
% 
% 
% % nesamone
% for i = 1:2*N
%     for j = 1:n
%         if X_cas(i,j) ~= 0 
%             E_total(i,j) = E_total(i,j) + E(j);
%         end
%     end
% end


E_total=zeros(2*N,3);
for i = 1:2*N
      for j =1:n
        if X_cas(i,j) == 1 
            E_total(i,1) =  E(j);
        elseif X_cas(i,j) == 2
                E_total(i,2) = E(j);
        elseif X_cas(i,j) == 3;
                    E_total(i,3) = E(j);
        end
      end

end


% initialise stress-induced death rate


% transition probability matrix excision
M_exc = zeros(N,N);

for i=1:n+1
   M_exc(i,1) = 1; 
end

% will use values of X_cas

for i = 1 : n*(n-1)
   M_exc(n+i+1,X_cas(n+i+1,1)+1) = 0.5;
   M_exc(n + i + 1, X_cas(n+i+1,2)+1) = 0.5;
end

n_2 = n + n*(n-1) + 2;


M_exc(n_2,5)=1/3;
M_exc(n_2,6)=1/3;
M_exc(n_2,8)=1/3;

M_exc(n_2+1,5)=1/3;
M_exc(n_2+1,6)=1/3;
M_exc(n_2+1,10)=1/3;

M_exc(n_2+2,6)=1/3;
M_exc(n_2+2,7)=1/3;
M_exc(n_2+2,8)=1/3;

M_exc(n_2+3,7)=1/3;
M_exc(n_2+3,8)=1/3;
M_exc(n_2+3,9)=1/3;

M_exc(n_2+4,5)=1/3;
M_exc(n_2+4,7)=1/3;
M_exc(n_2+4,10)=1/3;

M_exc(n_2+5,7)=1/3;
M_exc(n_2+5,9)=1/3;
M_exc(n_2+5,10)=1/3;


% matrix for re-integration

M_int = zeros(N,N);

for i=1:n+1
   M_int(i,i) = 1; 
end

   M_int(n+1+1,5) = 0.5;
   M_int(n + 1 + 1, 7) = 0.5;
   M_int(n+1+2,6) = 0.5;
   M_int(n + 1 + 2, 9) = 0.5;
   M_int(n+1+3,5) = 0.5;
   M_int(n + 1 + 3, 7) = 0.5;
   M_int(n+1+4,8) = 0.5;
   M_int(n + 1 + 4, 10) = 0.5;
   M_int(n+1+5,6) = 0.5;
   M_int(n + 1 + 5, 9) = 0.5;
   M_int(n+1+6,8) = 0.5;
   M_int(n + 1 + 6, 10) = 0.5;
   
   M_int(n + 1 + 7, 11) = 1/3;
   M_int(n + 1 + 7, 13) = 1/3;
   M_int(n + 1 + 7, 15) = 1/3;
   
   M_int(n + 1 + 8, 12) = 1/3;  S = [0;0;0];

   M_int(n + 1 + 8, 13) = 1/3;
   M_int(n + 1 + 8, 15) = 1/3;
   
   M_int(n + 1 + 9, 11) = 1/3;
   M_int(n + 1 + 9, 13) = 1/3;
   M_int(n + 1 + 9, 16) = 1/3;
  
   M_int(n + 1 + 10, 11) = 1/3;
   M_int(n + 1 + 10, 14) = 1/3;
   M_int(n + 1 + 10, 16) = 1/3;
   
   M_int(n + 1 + 11, 12) = 1/3;
   M_int(n + 1 + 11, 14) = 1/3;
   M_int(n + 1 + 11, 15) = 1/3;
   
   M_int(n + 1 + 12, 12) = 1/3;
   M_int(n + 1 + 12, 14) = 1/3;
   M_int(n + 1 + 12, 16) = 1/3;

%n_2 = n + n*(n-1) + 2;



%% solve ODEs




%initial condition
init_cond = zeros(N*2,1);
init_cond(11) = 10^6;






while time_step<10000
        
%     % lambda depends if the stressor is on or off
%         if S(i) == 0
%             lam = M(1,2);
%         else
%             lam = M(2,1); 
%         end
%         r = rand(1);
%         time = 1/lam * log(1/r); % different M(1,2), M(2,1) jei ijungtas ir jei isjungtas
        
      
% choose lambdas for each stressor depending on the states
for i = 1:3
     if S(i) == 0
            lam(1) = M(1,2);
            lam(2) = M(1,2);
            lam(3) = M(1,2);
     else
            lam(1) = M(2,1); 
            lam(2) = M(2,1);
            lam(3) = M(2,1);
     end

end       
       
% calculate the probability that one of three events happen
lam_prop(1) = lam(1)/sum(lam);
lam_prop(2) = lam(2)/sum(lam);
lam_prop(3) = lam(3)/sum(lam);

S_num = find( mnrnd(1,lam_prop));
        
    r = rand(1);
    time = 1/lam_prop(S_num) * log(1/r); % draw the next time step using expenential distribution
    
    
   S(S_num) =  mod(S(S_num)+1,2); % update the change in stress
             
        
    [t,y] = ode45(@(t,x) deriv(t,x,N,K,n_0,n_i,n_s,S,mu,ro,M_exc,M_int,theta,n,gamma,E_total),[time_step,time_step+time],init_cond);

    count = count+1;

    init_cond = y(size(t,1),:);
    time_step = time_step+time;
    
    t_final = [t_final;t];
    final = [final ; y];
    
    
  %  plot(t,y(:,12))
end

figure(1)
clf
nGenTypes = 16;

%% Plot results
% 
% figure(1)
% for g = 1:8
%     % Plot each genotype
%     subplot(8,1,g+1);
%     plot(t_final,final(:,g),'LineWidth',2,'LineStyle','-')
%     hold on
%     plot(t_final,final(:,nGenTypes+g),'LineWidth',2,'LineStyle','--')
%     hold off
%     legend('Xg','Yg')
%     title(['Genotype ' num2str(X_cas(g,:))])
% %     axis([0,T,0,K]);
% end
% shg

%plot(final)



figure(2)
% clf
% subplot(2,1,1)
% plot(t_final,SMat(:,1),'LineWidth',2,'LineStyle','-') 
% hold on
% plot(tVec,SMat(:,2),'LineWidth',2,'LineStyle','-')
% hold off2
% title('Stressors')
% legend('Stressor 1 (1=on, 0=off)','Stressor 2 (1=on, 0=off)')

subplot(2,1,2);
% Compute the total number of bacteria with gene i in first position
totPopwithGene1InFirstPos = final(:,2) + final(:,5) + final(:,6) + final(:,11) + final(:,12) +...
    + final(:,(2+nGenTypes)) + final(:,(5+nGenTypes))+ final(:,(6+nGenTypes))+final(:,(6+nGenTypes))+...
    +final(:,(11+nGenTypes))+final(:,(12+nGenTypes)); 
totPopwithGene2InFirstPos = final(:,3) + final(:,7)+ final(:,8)+ final(:,13)+ final(:,14)+...
    + final(:,(3+nGenTypes)) + final(:,(7+nGenTypes))+ final(:,(8+nGenTypes)) + final(:,(13+nGenTypes))+ final(:,(14+nGenTypes));
fracWithFunctIntegrase = sum(final(:,1:final),2)./sum(final,2);

%totPopwithGene3InFirstPos = final(:,3) + final(:,7)+ final(:,8)+ final(:,13)+ final(:,14)+...
%    + final(:,(3+nGenTypes)) + final(:,(7+nGenTypes))+ final(:,(8+nGenTypes)) + final(:,(13+nGenTypes))+ final(:,(14+nGenTypes));
%fracWithFunctIntegrase = sum(final(:,1:final),2)./sum(final,2);


% Plot it
plot(t_final,totPopwithGene1InFirstPos,'LineWidth',2,'LineStyle','-')
hold on
plot(t_final,totPopwithGene2InFirstPos,'LineWidth',2,'LineStyle',':')
% title(['Genotype ' num2str(genTypeMatrix(g,:))])

% yyaxis left
% %axis([0,T,0,K]);
% yyaxis right
% plot(t_final,fracWithFunctIntegrase,'LineWidth',2,'LineStyle','--')
% hold off
% legend('Gene 1 in First Position','Gene 2 in First Position', 'Proportion of total population with a functional Integrase')
% shg



    

