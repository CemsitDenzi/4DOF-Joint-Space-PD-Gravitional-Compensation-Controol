clc; clear; close all; 


t = 0:0.0001:10;
dt = 0.0001;
len = length(t);


q_ref = [0.2 * sin(0.2*t); 0.5 * sin(0.7*t); 0.6 * cos(0.4*t); 0.8 *cos(0.1*t)];
qd_ref = [0.2 * 0.2 * cos(0.2 * t);
          0.5 * 0.7 * cos(0.7 * t);
          -0.6 * 0.4 * sin(0.4 * t);
          -0.8 * 0.1 *sin(0.1 *t )];
qdd_ref = [-0.2 * 0.2 * 0.2 * sin(0.2 * t);
           -0.5 * 0.7 * 0.7 * sin(0.7 * t);
           -0.6 * 0.4 * 0.4 * cos(0.4 * t);
           -0.8 * 0.1 * 0.1 * cos(0.1 * t)];


torque = zeros(4, len); 

G = get_G(); 

%İnverse Dynamics
for i = 1:len

    q_i = q_ref(:, i);
    qd_i = qd_ref(:, i);
    
    D = get_D(q_i(1), q_i(2)); 
    C = get_C(q_i(1), q_i(2), qd_i(1), qd_i(2));
    
    torque(:, i) = D * qdd_ref(:, i) + C * qd_i + G;
end

q_sim = zeros(4, len);
qd_sim = zeros(4, len);

q_sim(:, 1) = q_ref(:, 1);
qd_sim(:, 1) = qd_ref(:, 1);

%Forward Dynamics
for i = 1:(len - 1) 

    D = get_D(q_sim(1, i), q_sim(2, i));
    C = get_C(q_sim(1, i), q_sim(2, i), qd_sim(1, i), qd_sim(2, i));
    G = get_G();

    qdd_res = inv(D) * (torque(:, i) - C * qd_sim(:, i) - G);

    qd_sim(:, i+1) = qd_sim(:, i) + qdd_res * dt;

    q_sim(:, i+1) = q_sim(:, i) + qd_sim(:, i+1) * dt;
end


figure;
for j=1:4
    subplot(4, 1, j);
    plot(t, q_ref(j, :), 'r--', 'LineWidth', 2);
    hold on;
    plot(t, q_sim(j, :), 'b-');
    grid on;
    title(['Eklem ', num2str(j), ' Pozisyonu']);
    legend('Referans (q_{ref})', 'Simülasyon (q_{sim})');
    xlabel('Zaman (s)');
    ylabel('Açı (rad)');
end