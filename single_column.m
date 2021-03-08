%% Parameters
N = 100;
tau = 0.001; %time constant
tau_ref = 0.003; %refractory period
tau_rec = 0.8
U = 0.5
J_ee = 6 %synaptic efficacy for exc to exc 
J_ei = -4 %synaptic efficacy for inh to exc 
J_ie = 0.5 %synaptic efficacy for exc to inh 
J_ii = -0.5 %synaptic efficacy for inh to inh 
e = -10 %background synaptic input
e_N = 10

%% Functions

    function dxdt = f(x)
             dxdt = (1-x/tau_rec)-U*x*E;
            
    end
    
    function dydt = f(y)
             dydt = (1-y/tau_rec)-U*y*I;
            
    end
    
    
    function dzdt = g(t,x)
        E = x(1);
        I = x(2);
        dzdt = [-E+(1-tau_ref*E)*((J_ee/N)*);
            -I+(1-tau_ref*I)];
    end