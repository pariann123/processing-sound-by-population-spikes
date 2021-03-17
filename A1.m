
function A1
N_E = 2;
N_I = 2;
P = 5;

global J_ee0 J_ei0 J_ii0 J_ie0 J_ee1 J_ee2 J_ie1 J_ie2;

J_ee0 = 6;
J_ei0 = -4;
J_ii0 = -0.5;
J_ie0 = 0.5;
J_ee1 = 0.045;
J_ee2 = 0.015;
J_ie1 = 0.0035;
J_ie2 = 0.0015;


lambdaX = 0.25;
alpha = 2;
delta_left = 5;
delta_right = 5;

vs= [];
for i=1:P
    vs0 =[rand(N_E,1);rand(N_I,1);rand(N_E,1);rand(N_I,1)]; % order is [E,I,x,y]
    vs = [vs, vs0];
end

rate_auditory(0,vs)

%% nested function
function out = rate_auditory(t,vs)

    E = vs(1:N_E,1:P);
    I = vs(N_E+1:N_E+N_I,1:P);
    x = vs(N_E+N_I+1: 2*N_E+N_I,1:P);
    y = vs(2*N_E+N_I+1:end,1:P);


    s=0; %sensory 
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
    
    test=[]
    for q = 3:P-2
        for R = -2:2
            var1 = j_ee(abs(R))/N_E;
            var2 = U.*x(1:N_E,q+R:q+R).*E(1:N_E,q+R:q+R); % check sum, it is summing by horizontal column    
            final = var1.*var2;
            test = [test,final];

        end
        sum1 = sum(test,2); %check if it should be sum horizontally.         
    end
     
    sum2 = J_ei0/N_I * sum(U.*y.*I);
    %s=0 BUT double check with markus
    
    out_E = max(0,sum1 + sum2 + ee_test); %relu
    
   
    
    


    sum_E = J_ee/N_E * sum(U.*x.*E) + (J_ei/N_I * (sum(U.*y.*I)));
    
    sum_I= J_ie/N_E * sum(E) + J_ii/N_I * sum(I);
    out_I = max(0,sum_I+ei_final); %relu

    dEdt = (-E + (1-tau_ref*E).*out_E)/tau;
    dIdt = (-I + (1-tau_ref*I).*out_I)/tau;
    dxdt = (1-x)/tau_rec - U*x.*E;
    dydt = (1-y)/tau_rec - U*y.*I;

    out = [dEdt; dIdt; dxdt; dydt]; % single column 
end


function out = j_ee(R)
        if R == 0
            out = 6;
        elseif R == 1
            out = 0.045;
        else
            out = 0.015;
        end
    
function out_h = heavf(x) %heaviside function    
    for i = 1:length(x)
        if x(i)>0
            out_h(i)=1;
        else
            out_h(i)=0;
        end
    end      
end
end
end



