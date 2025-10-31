clc;clear;

dt = 0.001;
simTime = 200 ; 

t = 0 : dt : simTime;
len = length(t);

q_ref = [0.2 * sin(0.2*t); 0.5 * sin(0.7*t); 0.6 * cos(0.4*t); 0.8 *cos(0.1*t)];
qd_ref = [0.2 * 0.2 * cos(0.2 * t);
          0.5 * 0.7 * cos(0.7 * t);
          -0.6 * 0.4 * sin(0.4 * t);
          -0.8 * 0.1 *sin(0.1 *t )];

q_res = zeros(4,len);
qd_res = zeros(4,len);

q_res(:,1) = q_ref(:,1);
qd_res(:,1) = qd_ref(:,1);


Kp = 100000;
Kd = 1000;
for(i = 1 : (len - 1))

    u = Kp*(q_ref(:,i) - q_res(:,i)) +  Kd * (qd_ref(:,i) - qd_res(:,i)) + get_G(); %Control signal [torque]
    D = get_D(q_res(1,i),q_res(2,i));
    C = get_C(q_res(1,i),q_res(2,i),qd_res(1,i),qd_res(2,i));
    G = get_G();

    qdd_res = inv(D) * (u - C * qd_res(:, i) - G);

    qd_res(:, i+1) = qd_res(:, i) + qdd_res * dt;

    q_res(:, i+1) = q_res(:, i) + qd_res(:, i+1) * dt;    


end

pos_err = q_ref - q_res;   

hold on;
plot(t, rad2deg(pos_err(1,:)))
plot(t, rad2deg(pos_err(2,:)))
plot(t, rad2deg(pos_err(3,:)))
plot(t, rad2deg(pos_err(4,:)))

grid on;
