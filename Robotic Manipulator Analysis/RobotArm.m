classdef RobotArm < handle

    properties
        l0
        l1
        l2
        l3
        dt
        Tf
        t
        Dt1
        Dt2
        q
        q_1
        J
        pE
        pE_1
    end

    methods
        function obj = RobotArm(l0, l1, l2, l3, dt, Tf)
            obj.l0 = l0;
            obj.l1 = l1;
            obj.l2 = l2;
            obj.l3 = l3;
            obj.dt = dt;
            obj.Tf = Tf;
            obj.t = 0:obj.dt:obj.Tf;
            duration = length(obj.t);
            obj.Dt1 = floor(duration/2);
            obj.Dt2 = duration - obj.Dt1;
            obj.pE = zeros(3, duration);
            obj.pE_1 = zeros(3, duration);
            obj.q = zeros(3, duration);
            obj.q_1 = zeros(3, duration);
            obj.J = zeros(6, 3, duration);
        end

        function generateTrajectory(o, PA, PB)
            % Positions
            Pc = (PB + PA) / 2;
            R = norm(Pc - PA);
            theta = linspace(0, pi, o.Dt2); % Only half the trajectory for the circle

            lin_pE = [
                linspace(PA(1), PB(1), o.Dt1);
                linspace(PA(2), PB(2), o.Dt1);
                linspace(PA(3), PB(3), o.Dt1);];

            cir_pE = [
                Pc(1) * ones(size(theta));
                Pc(2) + R*sin(theta);
                Pc(3) + R*cos(theta);];

            o.pE = [lin_pE, cir_pE];

            % Velocities
            lin_pE_1 = zeros(3, o.Dt1);
            for i = 1:3
                lin_pE_1(i, :) = (PB(i) - PA(i)) / (o.Dt1 * o.dt) * ones(1, o.Dt1);
            end
            omega = pi / (o.Dt2 * o.dt);
            cir_pE_1 = [
                zeros(size(theta));
                omega*R*cos(theta);
                -omega*R*sin(theta);];
            o.pE_1 = [lin_pE_1, cir_pE_1];
        end

        function calcJacobian(o)
            % Calculate Jacobian Matrix
            for k = 1:length(o.t)
                o.J(:,:,k) = [
                    -o.l1*sin(o.q(1, k)) + cos(o.q(1, k))*(o.l3*cos(o.q(2, k)+o.q(3, k)) - o.l2*sin(o.q(2, k))) , -sin(o.q(1, k))*(o.l3*sin(o.q(2, k)+o.q(3, k)) + o.l2*cos(o.q(2, k))) , -o.l3*sin(o.q(1, k))*sin(o.q(2, k)+o.q(3, k));
                    +o.l1*cos(o.q(1, k)) + sin(o.q(1, k))*(o.l3*cos(o.q(2, k)+o.q(3, k)) - o.l2*sin(o.q(2, k))) , +cos(o.q(1, k))*(o.l3*sin(o.q(2, k)+o.q(3, k)) + o.l2*cos(o.q(2, k))) , +o.l3*cos(o.q(1, k))*sin(o.q(2, k)+o.q(3, k));
                    0 , o.l3.*cos(o.q(2, k) + o.q(3, k)) - o.l2.*sin(o.q(2, k)) , o.l3.*cos(o.q(2, k)+o.q(3, k));
                    0 , -cos(o.q(1, k)) , -cos(o.q(1, k));
                    0 , +sin(o.q(1, k)) , +sin(o.q(1, k));
                    1 , 0 , 0
                    ];
            end
        end

        function performSimpleInverseKinematics(o)
            px = o.pE(1,:);
            py = o.pE(2,:);
            pz = o.pE(3,:);

            for k = 1:length(o.t)
                r = sqrt(px(k)^2 + py(k)^2);
                phi = atan2(py(k)/r, px(k)/r);

                o.q(1, k) = phi + atan2(sqrt(px(k)^2 + py(k)^2 - o.l1^2) / r, +o.l1 / r); % Positive solution for q(1)

                const1 = cos(o.q(1, k)) * py(k) - sin(o.q(1, k)) * px(k);
                const2 = pz(k) - o.l0;

                sigma = (const2^2 + const1^2 -o.l2^2 - o.l3^2) / (2 * o.l2 * o.l3);

                o.q(3, k) = atan2(sigma, +sqrt(1 - sigma^2)); % Positive solution for q(3)


                % Calculate q(2) based on q(1) and q(3)
                t1 = o.l2 + o.l3 * sin(o.q(3, k));
                t2 = o.l3 * cos(o.q(3, k));
                o.q(2, k) = atan2( ...
                    (t1*const1 + t2*const2) / (t1^2 + t2^2), ...
                    (t1*const2 - t2*const1) / (t1^2 + t2^2) ...
                    );

                o.J(:,:,k) = [
                    -o.l1*sin(o.q(1, k)) + cos(o.q(1, k))*(o.l3*cos(o.q(2, k)+o.q(3, k)) - o.l2*sin(o.q(2, k))) , -sin(o.q(1, k))*(o.l3*sin(o.q(2, k)+o.q(3, k)) + o.l2*cos(o.q(2, k))) , -o.l3*sin(o.q(1, k))*sin(o.q(2, k)+o.q(3, k));
                    +o.l1*cos(o.q(1, k)) + sin(o.q(1, k))*(o.l3*cos(o.q(2, k)+o.q(3, k)) - o.l2*sin(o.q(2, k))) , +cos(o.q(1, k))*(o.l3*sin(o.q(2, k)+o.q(3, k)) + o.l2*cos(o.q(2, k))) , +o.l3*cos(o.q(1, k))*sin(o.q(2, k)+o.q(3, k));
                    0 , o.l3.*cos(o.q(2, k) + o.q(3, k)) - o.l2.*sin(o.q(2, k)) , o.l3.*cos(o.q(2, k)+o.q(3, k));
                    0 , -cos(o.q(1, k)) , -cos(o.q(1, k));
                    0 , +sin(o.q(1, k)) , +sin(o.q(1, k));
                    1 , 0 , 0
                    ];

                % Find joint velocities
                o.q_1(:, k) = o.J(1:3,:,k) \ o.pE_1(:, k); % A\ => inv(A)
                % zero out the values less than 1e-10
                o.q_1(:, k) = o.q_1(:, k) .* (abs(o.q_1(:, k)) > 1e-10);
                o.q(:, k) = o.q(:, k) .* (abs(o.q(:, k)) > 1e-10);
            end
        end

        function [q, q_1] = calcJointState(o, pE, pE_1)
            q = zeros(3, 1);
            px = pE(1);
            py = pE(2);
            pz = pE(3);

            r = sqrt(px^2 + py^2);
            phi = atan2(py/r, px/r);

            q(1) = phi + atan2(sqrt(px^2 + py^2 - o.l1^2) / r, +o.l1 / r); % Positive solution for q(1)

            const1 = cos(o.q(1)) * py - sin(o.q(1)) * px;
            const2 = pz - o.l0;

            sigma = (const2^2 + const1^2 -o.l2^2 - o.l3^2) / (2 * o.l2 * o.l3);

            q(3) = atan2(sigma, +sqrt(1 - sigma^2)); % Positive solution for q(3)


            % Calculate q(2) based on q(1) and q(3)
            t1 = o.l2 + o.l3 * sin(q(3));
            t2 = o.l3 * cos(q(3));
            q(2) = atan2( ...
                (t1*const1 + t2*const2) / (t1^2 + t2^2), ...
                (t1*const2 - t2*const1) / (t1^2 + t2^2) ...
                );

            % Calculate Jacobian Matrix
            J_slice = [
                -o.l1*sin(q(1)) + cos(q(1))*(o.l3*cos(q(2)+q(3)) - o.l2*sin(q(2))) , -sin(q(1))*(o.l3*sin(q(2)+q(3)) + o.l2*cos(q(2))) , -o.l3*sin(q(1))*sin(q(2)+q(3));
                +o.l1*cos(q(1)) + sin(q(1))*(o.l3*cos(q(2)+q(3)) - o.l2*sin(q(2))) , +cos(q(1))*(o.l3*sin(q(2)+q(3)) + o.l2*cos(q(2))) , +o.l3*cos(q(1))*sin(q(2)+q(3));
                0 , o.l3.*cos(q(2) + q(3)) - o.l2.*sin(q(2)) , o.l3.*cos(q(2)+q(3));
                0 , -cos(q(1)) , -cos(q(1));
                0 , +sin(q(1)) , +sin(q(1));
                1 , 0 , 0
                ];

            % Find joint velocities
            q_1 = J_slice(1:3,:) \ pE_1; % A\ => inv(A) *
        end

        function [qi, qi_1, qf, qf_1] = calcIandFJointState(o, pi, pi_1, pf, pf_1)
            % Calculate initial and final joint state according to pi, pf.
            [qi, qi_1] = o.calcJointState(pi, pi_1);
            [qf, qf_1] = o.calcJointState(pf, pf_1);
        end

        function [qi, qf] = calcIandFJointAngles(o, pi, pf)
            % Calculate initial and final joint state according to pi, pf.
            [qi, ~] = o.calcJointState(pi, [0;0;0]);
            [qf, ~] = o.calcJointState(pf, [0;0;0]);
        end

        function [q, q_1] = calcJointTraj(o, qi, qi_1, qf, qf_1, ti, tf)
            % solution of A * o.a = bc for a => a = A^-1 * bc
            % Calculate coeficients
            duration = tf - ti;
            a(:, 1) = qi;
            a(:, 2) = qi_1;
            a(:, 3) = +(3/duration^2)*(qf-qi) - (2/duration)*qi_1 - qf_1/duration;
            a(:, 4) = -(2/duration^3)*(qf-qi) + (qi_1 + qf_1)/duration^2;

            % Calculate Joint State Values
            k = 0:o.dt:duration;
            A = [ones(size(k)); k; k.^2; k.^3];
            A_1 = [zeros(size(k)); ones(size(k)); 2*k; 3*k.^2];
            q = a * A;
            q_1 = a * A_1;
        end

        function smoothJointAnglesInverseKinematics(o)
            % For Linear trajectory i.e. PA to PB
            Dt1Seconds = o.Tf/2;
            Dt2Seconds = o.Tf - Dt1Seconds;
            [qiLin, qi_1Lin, qfLin, qf_1Lin] = o.calcIandFJointState(o.pE(:, 1), [0;0;0], o.pE(:, o.Dt1), [0;0;0]);
            [qLin, q_1Lin] = o.calcJointTraj(qiLin, qi_1Lin, qfLin, qf_1Lin, 0, Dt1Seconds);

            % For Circular Trajectory i.e. PB to PA split in 2 to do circle
            NSplits = 8;
            splitDuration = Dt2Seconds/NSplits - o.dt;
            qS_array = [];

            % Calculating q
            for iSplit = 0:(NSplits-1)
                tfIdx = floor(o.Dt1 + (iSplit+1) * o.Dt2 / NSplits);

                if (iSplit == 0)
                    [qi, qf] = o.calcIandFJointAngles(o.getpE(qLin(:, end)), o.pE(:, tfIdx));
                    qS_array = [qS_array, qi, qf];
                else
                    [~, qf] = o.calcIandFJointAngles(o.getpE(qS_array(:, end)), o.pE(:, tfIdx));
                    qS_array = [qS_array, qf];
                end
            end


            % Calculating q_1
            % Where: qf_1 = (sgn(vk) == sgn(vk_)) * (vk+v_k+1)/2  DIDN'T
            % MAKE SUCH A DIFFERENCE IN THE END...
            qS_1_array = zeros(3,1);
            qS_vel = diff(qS_array, 1, 2);
            for k = 1:NSplits-1
                qS_1_array = [qS_1_array, 0.5 * (sign(qS_vel(:, k)) == sign(qS_vel(:, k+1))) .* (qS_vel(:, k) + qS_vel(:, k+1))];
            end
            qS_1_array = [qS_1_array, zeros(3,1)];

            qCir = [];
            q_1Cir = [];
            for iSplit = 1:NSplits
                [qS, q_1S] = o.calcJointTraj(qS_array(:, iSplit), qS_1_array(:, iSplit), qS_array(:, iSplit+1), qS_1_array(:, iSplit+1), 0, splitDuration);
                qCir = [qCir, qS];
                q_1Cir = [q_1Cir, q_1S];
            end


            o.q = [qLin, qCir];
            o.q_1 = [q_1Lin, q_1Cir];
            % zero out the values less than 1e-10
            o.q_1 = o.q_1 .* (abs(o.q_1) > 1e-10);
            o.q = o.q .* (abs(o.q) > 1e-10);

            o.calcJacobian;
        end

        function plotJointAnglesInfo(o)
            figure;
            subplot(3,1,1);
            plot(o.t, o.q(1, :));
            title('Joint Angle q(1)');
            xlabel('Time (s)');
            ylabel('Angle (rad)');

            subplot(3,1,2);
            plot(o.t, o.q(2, :));
            title('Joint Angle q(2)');
            xlabel('Time (s)');
            ylabel('Angle (rad)');

            subplot(3,1,3);
            plot(o.t, o.q(3, :));
            title('Joint Angle q(3)');
            xlabel('Time (s)');
            ylabel('Angle (rad)');

            figure;
            subplot(3,1,1);
            plot(o.t, o.q_1(1,:));
            title('Joint Vel q(1)');
            xlabel('Time (s)');
            ylabel('Angle (rad)');

            subplot(3,1,2);
            plot(o.t, o.q_1(2,:));
            title('Joint Vel q(2)');
            xlabel('Time (s)');
            ylabel('Angle (rad)');

            subplot(3,1,3);
            plot(o.t, o.q_1(3,:));
            title('Joint Vel q(3)');
            xlabel('Time (s)');
            ylabel('Angle (rad)');
        end

        function pE = getpE(o, q_k)
            pE = [o.l1*cos(q_k(1)) + o.l3*sin(q_k(1))*cos(q_k(2)+q_k(3)) - o.l2*sin(q_k(1))*sin(q_k(2));
                o.l1*sin(q_k(1)) - o.l3*cos(q_k(1))*cos(q_k(2)+q_k(3)) + o.l2*cos(q_k(1))*sin(q_k(2));
                o.l0 + o.l2*cos(q_k(2)) + o.l3*sin(q_k(2)+q_k(3))];
        end

        function animateRobot(o, ax, showVelocities)
            px = o.pE(1,:);
            py = o.pE(2,:);
            pz = o.pE(3,:);

            pauseSeconds = 0.1;
            dtk = floor((pauseSeconds/o.Tf)*length(o.t)); %% plot robot position every dtk samples, to animate its motion

            p0 = [0;0;o.l0];
            plot3(ax, [0 p0(1)], [0 p0(2)], [0 p0(3)], 'y', 'MarkerSize', 10, 'MarkerFaceColor', 'k'); % Base
            hold (ax, "on");
            axis (ax, "equal");

            colors = ['r', 'g', 'b', 'm'];

            for k = 1:dtk:length(o.t)
                pause(pauseSeconds);	%% pause motion to view successive robot configurations
                % Calculate positions based on joint angles
                p1_k = [o.l1*cos(o.q(1, k)); o.l1*sin(o.q(1, k)); o.l0];
                p2_k = [o.l1*cos(o.q(1, k)) - o.l2*sin(o.q(1, k))*sin(o.q(2, k));
                    o.l1*sin(o.q(1, k)) + o.l2*cos(o.q(1, k))*sin(o.q(2, k));
                    o.l2*cos(o.q(2, k)) + o.l0];

                pE_k = o.getpE(o.q(:, k));

                % Plotting the robot's segments
                plot3(ax, [p0(1) p1_k(1)], [p0(2) p1_k(2)], [p0(3) p1_k(3)], colors(1), 'LineWidth', 2); % Base to Joint 1
                plot3(ax, [p1_k(1) p2_k(1)], [p1_k(2) p2_k(2)], [p1_k(3) p2_k(3)], colors(2), 'LineWidth', 2); % Joint 1 to Joint 2
                plot3(ax, [p2_k(1) pE_k(1)], [p2_k(2) pE_k(2)], [p2_k(3) pE_k(3)], colors(3), 'LineWidth', 2); % Joint 2 to End Effector

                plot3(ax, px, py, pz)

                % Plot end effector position
                plot3(ax, pE_k(1), pE_k(2), pE_k(3), 'm*', 'MarkerSize', 5);




                Jp_k = o.J(1:3, :, k);
                pE_1_k = Jp_k * o.q_1(:, k);

                if showVelocities
                    quiver3(ax, pE_k(1), pE_k(2), pE_k(3), pE_1_k(1), pE_1_k(2), pE_1_k(3));
                end
            end
            hold (ax, "off");
        end
    end
end
