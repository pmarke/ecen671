Ts = 0.01;
system = systemDynamics(Ts);  
rls = rlsDynamics(Ts);
y_history = [];
u_history = [];
t_history = [];
a1_history = [];
a0_history = [];
b0_history = [];
x_rls = [0; 0; 0];

% main simulation loop
t = 0;  % time starts at t=0
while t < 100  
    y = system.output;
    u = sin(0.01*t)+sin(0.1*t)+sin(t)+sin(10*t);
    system.propagateDynamics(u);  % Propagate the dynamics
    
    % update rls algorithm
    rls.update(y, u);
    x_rls = rls.state;
    
    % update histories
    y_history = [y_history; y];    
    u_history = [u_history; u];
    t_history = [t_history; t];
    a1_history = [a1_history; x_rls(1)];
    a0_history = [a0_history; x_rls(2)];
    b0_history = [b0_history; x_rls(3)];
    
    % advance time by Ts
    t = t + Ts;     
end
[a1_ls, a0_ls, b0_ls] = batch_ls(y_history, u_history, Ts);

% create plot
figure(1), clf
    subplot(311)
    plot(t_history, a1_history,'b')
    hold on
    plot(t_history, a1_ls*ones(size(t_history)),'g')
    ylabel('a_1')

    subplot(312)
    plot(t_history, a0_history,'b')
    hold on
    plot(t_history, a0_ls*ones(size(t_history)),'g')
    ylabel('a_0')

    subplot(313)
    plot(t_history, b0_history,'b')
    hold on
    plot(t_history, b0_ls*ones(size(t_history)),'g')
    ylabel('b_0')

