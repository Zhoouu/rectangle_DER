%% DDG tutorial, 3d_rod, Case 2: Bifurcation of a pre-buckled ribbon
% Weicheng Huang, weicheng.huang@ncl.ac.uk
% Dezhong Tong, dezhong@umich.edu
% Zhuonan Hao, znhao@g.ucla.edu

clear all;
close all;
clc;

global global_x0;
% Discrete ribbon simulation
fprintf('Rectangle Rotating Rod Under Gravity Contact\n');

% input nodes
node = importdata('inputfile/two_rod_nodes.txt');

% input stretching element
edge = importdata('inputfile/two_rod_edge.txt');

% input bending element
bend = importdata('inputfile/two_rod_bends.txt');

% Numerical parameter
simParams = defSimParams();

% Build beam struct
rodParams = defRodParams(node, edge, bend, simParams);
global_x0=rodParams.x0;
% Build stretching element
sElement = InitialStretchingElement(rodParams, edge);

% Build bending element
bElement = InitialBendingElement(rodParams, bend, sElement);

% Build boundary condition
consParams = defConsParams(rodParams);

% Current time
ctime = 0.0;

% --- 1. 定义边的颜色 ---
numEdges = size(edge, 1); % 获取边的总数

% 首先，创建一个元胞数组，将所有边的默认颜色设置为黑色 ('k')
edgeColors = repmat({'k'}, numEdges, 1);

middle_edges_to_color_red=99:123;
edgeColors(middle_edges_to_color_red) = {'r'}; % 将指定的多条边设置为红色

% Get constrained dof and unconstrained dof
rodParams.xCons  = rodParams.x(consParams.consInd);
rodParams.xUncons = rodParams.x(consParams.unconsInd);

%% 设置第一步压缩屈曲参数
totalCompress = 0.0;
buckling_time_start = 0.01; % Time to start compression
buckling_time_end = 2;   % Time to end compression/buckling phase
targetCompress = 0.05; % Target compression amount (adjust as needed)
compression_speed = targetCompress / (buckling_time_end - buckling_time_start);

%% 设置第二步旋转参数
fixed_end_end_position=289:300;
rotating_implement = [fixed_end_end_position'];
separate_dist=global_x0(3);
rotation_start_time = buckling_time_end; % Start immediately after buckling
rotation_end_time = 4.5;
total_rotation_angle = -pi/2*2.6; % Target rotation angle in radians (e.g., -90 degrees) Adjust as needed
rotation_duration = rotation_end_time - rotation_start_time;
angular_velocity = total_rotation_angle / rotation_duration;
current_angle = 0; % Start rotation from the buckled state's segment angle
% Open file for writing
fileID = fopen('data.txt', 'w'); 

%% 设置第三步压缩屈曲参数
totalCompress_2 = targetCompress;
buckling_time_start_2 = rotation_end_time; % Time to start compression
buckling_time_end_2 = simParams.totalTime;   % Time to end compression/buckling phase
targetCompress_2 = 0.2; % Target compression amount (adjust as needed)
compression_speed_2 = (targetCompress_2 - totalCompress_2) / (buckling_time_end_2 - buckling_time_start_2);

% --- 添加初始化代码 ---
all_iters = zeros(simParams.Nsteps, 1);  % 存储每一步的迭代次数
all_errors = zeros(simParams.Nsteps, 1); % 存储每一步的最终误差
max_iter_overall = 0;                    % 记录整个过程中的最大迭代次数
max_error_overall = 0;                   % 记录整个过程中的最大误差
max_error_step = 0;                      % 记录最大误差发生在哪一步
% --- 结束初始化代码 ---

% Shear a ribbon to induce bifurcation
totalShear = 0.0;
inputShear = 0.0;
%% --- GIF Initialization ---
gifFilename = ['buckling_rotation_buckling.gif']; % Name of the output GIF file
frameDelay = 0.1;      % Delay between frames in seconds (adjust as needed)
isFirstFrame = true;   % Flag to handle the first frame write differently
fprintf('Initializing GIF generation. Output file: %s\n', gifFilename);
% Create a figure window specifically for the animation BEFORE the loop
hFig_gif = figure; 
set(hFig_gif, 'Name', 'two rod rotate'); % Optional: Set figure name

%% Simulation loop
buckling_steps = round(buckling_time_end / simParams.dt);

for timeStep=1:simParams.Nsteps
    
    fprintf('t=%f\n', ctime);
    
    if (ctime <= buckling_time_start)
        rodParams.g = [0.0;0; -1];
        rodParams.ifStatic = 0;  % 为了添加外部力项
    end

    % Compress the ribbon
    if (ctime > buckling_time_start && totalCompress <=targetCompress)
        for ii=1:2
            rodParams.x(3*(ii-1)+1) = rodParams.x(3*(ii-1)+1) + simParams.dt * compression_speed;
            rodParams.x(3*(ii-1)+7) = rodParams.x(3*(ii-1)+7) + simParams.dt * compression_speed;
        end
        % rodParams.x(3*rodParams.nv-2) = rodParams.x(3*rodParams.nv-2) - simParams.dt * 0.1;
        % rodParams.x(3*rodParams.nv-5) = rodParams.x(3*rodParams.nv-5) - simParams.dt * 0.1;
        totalCompress = totalCompress + simParams.dt * compression_speed;
    end

    if(ctime > buckling_time_start+simParams.dt)
        rodParams.g=[0.0;0;-1];
        rodParams.ifStatic = 0; % 为了添加外部力项
    end
    % 
    % % Shear the ribbon
    % if (ctime > 5.0 && totalShear <= inputShear)
    % 
    %     rodParams.x(2) = rodParams.x(2) + simParams.dt * 0.01;
    %     rodParams.x(5) = rodParams.x(5) + simParams.dt * 0.01;
    %     rodParams.x(3*rodParams.nv-1) = rodParams.x(3*rodParams.nv-1) - simParams.dt * 0.01;
    %     rodParams.x(3*rodParams.nv-4) = rodParams.x(3*rodParams.nv-4) - simParams.dt * 0.01;
    % 
    %     totalShear = totalShear + 2 * simParams.dt * 0.01;
    % end
    
    if (ctime > rotation_start_time&&ctime<buckling_time_start_2)
        % rodParams.x(3*rodParams.nv+1)=current_angle;
        current_angle=angular_velocity * (ctime - rotation_start_time);
        % rodParams.x(2) = separate_dist*sin(current_angle);
        % rodParams.x(3) = separate_dist*cos(current_angle);
        % rodParams.x(5) = separate_dist*sin(current_angle);
        % rodParams.x(6) = separate_dist*cos(current_angle);
        R=[cos(current_angle) -sin(current_angle);sin(current_angle) cos(current_angle)];
        P=Position_end_buckle_end*R;
        rodParams.x(290)=P(1,1);
        rodParams.x(291)=P(1,2);
        rodParams.x(293)=P(2,1);
        rodParams.x(294)=P(2,2);
        rodParams.x(296)=P(3,1);
        rodParams.x(297)=P(3,2);
        rodParams.x(299)=P(4,1);
        rodParams.x(300)=P(4,2);
        % rodParams.x(1178)=P(1,1);
        % rodParams.x(1179)=P(1,2);
        % rodParams.x(1181)=P(2,1);
        % rodParams.x(1182)=P(2,2);
        % rodParams.x(1184)=P(3,1);
        % rodParams.x(1185)=P(3,2);
        % rodParams.x(1187)=P(4,1);
        % rodParams.x(1188)=P(4,2);
        % rodParams.x(1190)=P(5,1);
        % rodParams.x(1191)=P(5,2);
        % rodParams.x(1193)=P(6,1);
        % rodParams.x(1194)=P(6,2);
        % rodParams.x(1196)=P(7,1);
        % rodParams.x(1197)=P(7,2);
        % rodParams.x(1199)=P(8,1);
        % rodParams.x(1200)=P(8,2);
        % rodParams.x(1154)=P(9,1);
        % rodParams.x(1155)=P(9,2);
        % rodParams.x(1157)=P(10,1);
        % rodParams.x(1158)=P(10,2);
        % rodParams.x(1160)=P(11,1);
        % rodParams.x(1161)=P(11,2);
        % rodParams.x(1163)=P(12,1);
        % rodParams.x(1164)=P(12,2);
        % rodParams.x(1166)=P(13,1);
        % rodParams.x(1167)=P(13,2);
        % rodParams.x(1169)=P(14,1);
        % rodParams.x(1170)=P(14,2);
        % rodParams.x(1172)=P(15,1);
        % rodParams.x(1173)=P(15,2);
        % rodParams.x(1175)=P(16,1);
        % rodParams.x(1176)=P(16,2);
    end

    % Compress the ribbon
    if (ctime >= buckling_time_start_2 && totalCompress_2 <=targetCompress_2)
        for ii=1:2
            rodParams.x(3*(ii-1)+1) = rodParams.x(3*(ii-1)+1) + simParams.dt * compression_speed_2;
            rodParams.x(3*(ii-1)+7) = rodParams.x(3*(ii-1)+7) + simParams.dt * compression_speed_2;
        end
        % rodParams.x(3*rodParams.nv-2) = rodParams.x(3*rodParams.nv-2) - simParams.dt * 0.1;
        % rodParams.x(3*rodParams.nv-5) = rodParams.x(3*rodParams.nv-5) - simParams.dt * 0.1;
        totalCompress_2 = totalCompress_2 + simParams.dt * compression_speed_2;
    end
    % Initial guess
    rodParams.x(consParams.unconsInd) = rodParams.x(consParams.unconsInd) + rodParams.u(consParams.unconsInd) * rodParams.dt;

    % Solver
    [xUncons, iter_this_step, normf_this_step] = objfun(rodParams, simParams, consParams, sElement, bElement, edge);
    
    all_iters(timeStep) = iter_this_step;
    all_errors(timeStep) = normf_this_step;
     % --- 添加记录和更新最大值的代码 ---
    all_iters(timeStep) = iter_this_step;
    all_errors(timeStep) = normf_this_step;

    if iter_this_step > max_iter_overall
        max_iter_overall = iter_this_step;
    end

    if normf_this_step > max_error_overall
        max_error_overall = normf_this_step;
        max_error_step = timeStep; % 记录发生最大误差的时间步
    end
    % --- 结束记录代码 ---

    % Update DOF
    rodParams.x(consParams.unconsInd) = xUncons;
    
    % Update velocity
    rodParams.u = (rodParams.x - rodParams.x0) / rodParams.dt;
    
    % Update x0
    rodParams.x0 = rodParams.x;
    
    % Update time
    ctime = ctime + simParams.dt;
    
    % Update stretching element 
    for i = 1:rodParams.ne
        sElement(i).nodePos_1 = getVertex(rodParams.x, sElement(i).nodeIndex(1));
        sElement(i).nodePos_2 = getVertex(rodParams.x, sElement(i).nodeIndex(2));
        edge_vec = sElement(i).nodePos_2 - sElement(i).nodePos_1;
        sElement(i).t = edge_vec / norm(edge_vec);
    end
    
    for i =1:rodParams.nb
        % Store old frames before updating
        bElement(i).t_1_old  = bElement(i).t_1;
        bElement(i).d_11_old = bElement(i).d_11;
        bElement(i).d_12_old = bElement(i).d_12;
        bElement(i).t_2_old  = bElement(i).t_2;
        bElement(i).d_21_old = bElement(i).d_21;
        bElement(i).d_22_old = bElement(i).d_22;
        bElement(i).refTwist_old = bElement(i).refTwist; % Store old refTwist

        % Update positions
        bElement(i).nodePos_1 = getVertex(rodParams.x, bElement(i).nodeIndex(1));
        bElement(i).nodePos_2 = getVertex(rodParams.x, bElement(i).nodeIndex(2));
        bElement(i).nodePos_3 = getVertex(rodParams.x, bElement(i).nodeIndex(3));

        % Update theta
        if (bElement(i).directSign_1 > 0)
            bElement(i).theta_1 = getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
        else
            bElement(i).theta_1 = - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
        end
        if (bElement(i).directSign_2 > 0)
            bElement(i).theta_2 = getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
        else
            bElement(i).theta_2 = - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
        end

        % Recompute frames and twist based on converged state
        t1_new = (bElement(i).nodePos_2 - bElement(i).nodePos_1);
        t2_new = (bElement(i).nodePos_3 - bElement(i).nodePos_2);
        norm_t1 = norm(t1_new); norm_t2 = norm(t2_new);
        if norm_t1 < 1e-12 || norm_t2 < 1e-12
             % Keep old frames if edge length is near zero (already warned in objfun)
             bElement(i).t_1 = bElement(i).t_1_old;
             bElement(i).t_2 = bElement(i).t_2_old;
             bElement(i).d_11 = bElement(i).d_11_old;
             bElement(i).d_12 = bElement(i).d_12_old;
             bElement(i).d_21 = bElement(i).d_21_old;
             bElement(i).d_22 = bElement(i).d_22_old;
             bElement(i).refTwist = bElement(i).refTwist_old;
        else
            bElement(i).t_1 = t1_new / norm_t1;
            bElement(i).t_2 = t2_new / norm_t2;

            % Update directors via parallel transport from OLD frames to NEW tangents
            bElement(i).d_11 = parallel_transport(bElement(i).d_11_old, bElement(i).t_1_old, bElement(i).t_1);
            bElement(i).d_11 = bElement(i).d_11 / norm(bElement(i).d_11);
            bElement(i).d_12 = cross(bElement(i).t_1, bElement(i).d_11);

            bElement(i).d_21 = parallel_transport(bElement(i).d_21_old, bElement(i).t_2_old, bElement(i).t_2);
            bElement(i).d_21 = bElement(i).d_21 / norm(bElement(i).d_21);
            bElement(i).d_22 = cross(bElement(i).t_2, bElement(i).d_21);

            % Update reference twist based on the change from OLD frames/twist
            u_temp = parallel_transport(bElement(i).d_11, bElement(i).t_1, bElement(i).t_2);
            u_temp = rotateAxisAngle(u_temp, bElement(i).t_2, bElement(i).refTwist_old);
            deltaAngle = signedAngle(u_temp, bElement(i).d_21, bElement(i).t_2);
            bElement(i).refTwist = bElement(i).refTwist_old + deltaAngle; % FINAL updated twist for this step
        end

        % Update material frames
        cs = cos(bElement(i).theta_1); ss = sin(bElement(i).theta_1);
        bElement(i).m_11 = cs * bElement(i).d_11 + ss * bElement(i).d_12;
        bElement(i).m_12 = -ss * bElement(i).d_11 + cs * bElement(i).d_12;

        cs = cos(bElement(i).theta_2); ss = sin(bElement(i).theta_2);
        bElement(i).m_21 = cs * bElement(i).d_21 + ss * bElement(i).d_22;
        bElement(i).m_22 = -ss * bElement(i).d_21 + cs * bElement(i).d_22;
    end
    % --- End Element Data Update ---

    if (ctime<=buckling_time_end)
          Position_end_buckle_end=[rodParams.x(290) rodParams.x(291);rodParams.x(293) rodParams.x(294);rodParams.x(296) rodParams.x(297);rodParams.x(299) rodParams.x(300)]; 
    end
    % % Update bending element 
    % for i =1:rodParams.nb
    %     bElement(i).nodePos_1 = getVertex(rodParams.x, bElement(i).nodeIndex(1));
    %     bElement(i).nodePos_2 = getVertex(rodParams.x, bElement(i).nodeIndex(2));
    %     bElement(i).nodePos_3 = getVertex(rodParams.x, bElement(i).nodeIndex(3));
    % 
    %     if (bElement(i).directSign_1 > 0)
    %         bElement(i).theta_1 =    getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
    %     else
    %         bElement(i).theta_1 =  - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(1));
    %     end
    % 
    %     if (bElement(i).directSign_2 > 0)
    %         bElement(i).theta_2 =    getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
    %     else
    %         bElement(i).theta_2 =  - getTheta(rodParams.x, rodParams.nv, bElement(i).edgeIndex(2));
    %     end
    % 
    %     bElement(i).t_1_old  = bElement(i).t_1;
    %     bElement(i).d_11_old = bElement(i).d_11;
    %     bElement(i).d_12_old = bElement(i).d_12; 
    % 
    %     bElement(i).t_2_old  = bElement(i).t_2;
    %     bElement(i).d_21_old = bElement(i).d_21;
    %     bElement(i).d_22_old = bElement(i).d_22;
    % 
    %     bElement(i).refTwist_old = bElement(i).refTwist;
    % 
    %     cs = cos( bElement(i).theta_1 );
    %     ss = sin( bElement(i).theta_1 );
    %     bElement(i).m_11 =   cs * bElement(i).d_11 + ss * bElement(i).d_12;
    %     bElement(i).m_12 = - ss * bElement(i).d_11 + cs * bElement(i).d_12;
    % 
    %     cs = cos( bElement(i).theta_2 );
    %     ss = sin( bElement(i).theta_2 );
    %     bElement(i).m_21 =   cs * bElement(i).d_21 + ss * bElement(i).d_22;
    %     bElement(i).m_22 = - ss * bElement(i).d_21 + cs * bElement(i).d_22;
    % end
    %% 开始绘制动图，注意前面还需要进行 Gif Initialization
    if (mod(timeStep-1, simParams.plotStep) == 0)
        % Make sure plotting happens in our dedicated figure
        figure(hFig_gif); 
        
        % Adjust plot limits if needed (can be done inside plotBeam or here)
        % plotRod(rodParams.x, rodParams.nv); % Call your plotting function
        plotRod(rodParams.x, rodParams.nv, edge, edgeColors); % Call your plotting function
        title(sprintf('Rotate State at t=%.3f s (Step %d)', ctime, timeStep)); % Add informative title
        
        % Force MATLAB to draw the plot NOW
        drawnow; 
        
        % Capture the frame from the dedicated figure
        frame = getframe(hFig_gif);
        img = frame2im(frame);
        [imind, cm] = rgb2ind(img, 256); % Convert to indexed image for GIF

        % Write frame to GIF
        if isFirstFrame
            imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', frameDelay);
            isFirstFrame = false; % Update flag
        else
            imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
        end
        fprintf('    Captured frame %d for GIF.\n', timeStep); % Progress indicator
    end  

    %% Plot figure
    if (mod(timeStep, simParams.plotStep) == 0)
        % plotRod(rodParams.x, rodParams.nv);
        plotRod(rodParams.x, rodParams.nv, edge, edgeColors); % Call your plotting function
    end
    
    % Cout data
    if (ctime > 5.0)
         for i = 1:rodParams.nb
            localNode = bElement(i).nodePos_2;
            m1 = (bElement(i).m_11 + bElement(i).m_21) / 2;
            outData = [localNode;m1];
            fprintf(fileID, '%.4f %.4f %.4f %.4f %.4f %.4f %.4f \n', [totalShear outData']);  % Custom formatting
         end
    end
end


fprintf('\n--- Simulation Summary ---\n');
fprintf('Total time steps: %d\n', simParams.Nsteps);
fprintf('Maximum Newton iterations in any time step: %d\n', max_iter_overall);
fprintf('Maximum final residual error in any time step: %.4e (occurred at step %d)\n', max_error_overall, max_error_step);
fprintf('Simulation finished. GIF saved as %s\n', gifFilename);

fclose(fileID);
