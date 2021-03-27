function A1
N_E = 2;
N_I = 2;
P = 5;
E_range = N_E;
I_range = N_E+N_I;
x_range = 2*N_E+N_I;
lambdaC = 0.25;
alpha = 2;
D = 5
D_left = 5;
D_right = 5;
J_ei0 = -4;
J_ii0 = -0.5;


z = 0; %sensory
U = 0.5;
tau = 0.001; % tau constant
tau_ref = 0.003;
tau_rec = 0.8;


%% background synaptic inputs
% excitatory
ee_test = [];
ei_test = [];

for i = 1:P
    e_temp = rand(N_E,1); % initial random sampling of N_E
    e_temp_min = min(e_temp);
    e_temp_max = max(e_temp);
    e_min = -10; % parameters
    e_max = 10; % parameters
    e = ((e_max-e_min)*(e_temp-e_temp_min)/(e_temp_max-e_temp_min))+e_min; % scaling to e_min + e_max
    ee_final = sort(e); % uniform distribution
    ee_test = [ee_test,ee_final];
    
    % inhibitory
    e_temp = rand(N_I,1); % initial random sampling of N_I
    e_temp_min = min(e_temp);
    e_temp_max = max(e_temp);
    e_min = -10;
    e_max = 10;
    e = ((e_max-e_min)*(e_temp-e_temp_min)/(e_temp_max-e_temp_min))+e_min; % scaling
    ei_final = sort(e);
    ei_test = [ei_test,ei_final];
end


vs= [];
for i=1:P
    vs0 =[rand(N_E,1);rand(N_I,1);rand(N_E,1);rand(N_I,1)]; % order is [E,I,x,y]
    vs = [vs; vs0];
end

rate_auditory(0,vs)
tspan = [0 0.1];
[tt,xx] = ode45(@rate_auditory,tspan,vs);


OE = xx(1:length(tt),1:E_range*P);
OI = xx(1:length(tt),E_range*P+1:I_range*P);
Ox = xx(1:length(tt),I_range*P+1: x_range*P);
Oy = xx(1:length(tt),x_range*P+1:end);

E_range = N_E;
I_range = N_E + N_I;
x_range = 2*N_E + N_I;

I = vs(E_range*P+1:I_range*P);
        I_mat = reshape(I, N_I, P);
        x = vs(I_range*P+1:x_range*P);
        x_mat = reshape(x, N_E,P);
        y = vs(x_range*P+1:end);
        y_mat = reshape(y, N_I,P);
        
% time, columns, neurons
mOE = zeros(length(tt),P);
mOI = zeros (length(tt),P);
mOx = zeros(length(tt),P);
mOy = zeros(length(tt),P);

for i=1:P
    mOE(:,i) = mean(OE(:,1:E_range),2);
    mOI(:,i) = mean(OI(E_range+1:I_range),2);
    mOx(:,i) = mean(Ox(I_range+1:x_range),2);
    mOy(:,i) = mean(Oy(x_range+1:end),2);
end

mOE = mean(OE,2);
mOI = mean(OI,2);
close all;
figure;
it = floor(length(tt)/2);
plot(mOE,mOI);
hold on
title('phase-plane diagram')

figure;
subplot(2,2,1); plot(tt,OE); title('Excitatory')
hold on; plot(tt,mOE,'linewidth',2);
subplot(2,2,2); plot(tt,OI); title('Inhibitory')
subplot(2,2,3); plot(tt,Ox); title('x')
subplot(2,2,4); plot(tt,Oy); title('y')

% %% nested function
    function out = rate_auditory(t,vs)
        
        % state variables in matrix and vector form
        
        E = vs(1:E_range*P);
        E_mat = reshape(E,N_E,P);
        I = vs(E_range*P+1:I_range*P);
        I_mat = reshape(I, N_I, P);
        x = vs(I_range*P+1:x_range*P);
        x_mat = reshape(x, N_E,P);
        y = vs(x_range*P+1:end);
        y_mat = reshape(y, N_I,P);
        
        if t > 0
            z = 4; %zeta
        end
        
        
        for q=1:P
            switch q
                case 1
                    R_range = 0:2;
                case 2
                    R_range = -1:2;
                case P-1
                    R_range = -2:1;
                case P
                    R_range = -2:0;
                otherwise
                    R_range = -2:2;
            end
            q_sumE=[];
            q_sumI = [];
            sum1_E=[];
            sum1_I= [];
            for R = R_range
                var1 = j_ee(abs(R))/N_E;
                var2 = sum(U*x_mat(:,q+R).*E_mat(:,q+R));
                final = var1.*var2;
                q_sumE = [q_sumE,final];
                sum_e = sum(E_mat(:,q+R)); %sum of E for the I rate
                sum_I = (j_ie(abs(R))/N_E) * sum_e;
                q_sumI = [q_sumI,sum_I];
            end
            sum1_E = [sum1_E,sum(q_sumE)];
            sum1_I = [sum1_I, sum(q_sumI)];
            
            %h = ones(P,1); % change h
%             h = 2*(abs(P- P/2))/P *ones(P,1); % middle neuron receives the strongest
            h = sptial(A)*ones(P,1)
            s = z*h;
            sum3_E(q) = sum(s);
        end
        
        sum2_E = (J_ei0/N_I) * sum(U.*y_mat.*I_mat);
        %s=0 BUT double check with markus
        
        
        out_E = max(0,sum1_E + sum2_E + ee_test + sum3_E); %relu
        
        
        sum_I = sum1_I +J_ii0/N_I * sum(I_mat);
        out_I = max(0,sum_I+ei_final); %relu
        
        dEdt = (-E_mat + (1-tau_ref*E_mat).*out_E)/tau;
        dIdt = (-I_mat + (1-tau_ref*I_mat).*out_I)/tau;
        dxdt = (1-x_mat)/tau_rec - U*x_mat.*E_mat;
        dydt = (1-y_mat)/tau_rec - U*y_mat.*I_mat;
        
        out= [dEdt(:);dIdt(:);dxdt(:);dydt(:)];
    end


    function out = j_ee(R)
        if R == 0
            out = 6;
        elseif R == 1
            out = 0.045;
        else
            out = 0.015;
        end
    end

    function out = j_ie(R)
        if R == 0
            out = 6;
        elseif R == 1
            out = 0.0035;
        else
            out = 0.0015;
        end
    end

    function out = spatial(A)
%         for i = 1:length(x)
            if A <= alpha
                peak_mag = lambdaC;
            else
                peak_mag = lambdaC + ((A-alpha)/D; %do we need D_left?
            end      
            
        out = A*exp(-abs(Q-M)/peak_mag;
    end
end