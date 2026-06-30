% ============================================================
% 静基座平台式惯导/里程计组合仿真单文件版本
% ============================================================

% 清空工作区、命令行和图窗，避免旧变量或旧图影响本次仿真。
clear; clc; close all;
% 初始化全局常量：地球半径、自转角速度、重力、角度单位、陀螺/加速度计单位等。
gvar;
% 这里只在主脚本中用到 arcdeg，用于把 30 deg 转成弧度。
global arcdeg
% 运行一个静基座场景：
%   sceneName = 'Static base'：结果打印和图标题使用的场景名称；
%   vnRef = [0;0;0]：参考速度为东、北、天三方向均为 0，即静止；
%   attb0 = [0;0;30]*arcdeg：车体系相对导航系存在 30 deg 航向角；
%   seed = 11：固定随机种子，保证加入的噪声每次一致。
initialAlignAtt = [0; 0; 0];  % [roll; pitch; yaw], rad; simulation route sets the heading internally.
result_static = platform_odo_real_imu_core('Beijing 4th-ring rectangle simulation', [0; 0; 0], initialAlignAtt, 11);
% 保存 result_static 结构体，里面包含误差曲线、滤波状态、传感器误差参数和均方误差。
save('platform_odo_fourth_ring_sim_results.mat', 'result_static');

% ============================================================
% 核心仿真函数：完成一次平台惯导 + 里程计组合导航仿真
% 输入：
%   sceneName：场景名称，仅用于显示；
%   vnRef：真实/参考速度，排列为 [vE; vN; vU]，即东、北、天；
%   attb0：载体姿态角，作为里程计坐标分解矩阵 Cbn0 的来源；
%   seed：随机数种子，用于复现实验噪声。
% 输出 result：保存误差、状态估计、噪声参数和统计指标。
% ============================================================
function result = platform_odo_real_imu_core(sceneName, vnRef, attb0, seed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 模块一：初始化模块
% 输入：sceneName、vnRef、attb0、seed
% 输出：地理参数、滤波器参数、惯性器件误差参数、惯导初始状态
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% （1）地理参数和仿真基础参数初始化
% 在函数内部再次初始化全局常量，保证单独调用该函数时也能正常运行。
    gvar;
% arcmin/dph/ug 等单位常量用于设置初始误差、传感器误差和绘图换算；Re 用于经纬度误差换算成米。
    global arcdeg arcmin dph dpsh ug ugpsHz Re g0
% 固定随机数序列，使 IMU 随机游走和里程计噪声在多次运行时完全可重复。
    rng(seed);
% nn 是每次惯导更新合并的 IMU 小采样数；ts 是 IMU 采样周期；nts 是一次导航解算更新周期。
% 这里 nn=2、ts=0.1 s，所以导航状态每 0.2 s 更新一次。
    ts = 0.1;
    nn = 1;
    nts = nn*ts;
% 真实 IMU 的 MAT 文件接口。每次导航更新读取 50 行、每行采样周期为 2 ms。
    imuCfg.useRealData = false;
    imuCfg.source = 'mat';
    imuCfg.file = 'dataCS1034-C09#自瞄准好 HSP-2026-6-5-13-38-43_part1.mat';
% 留空时，程序会自动选择 MAT 文件中第一个至少包含六列的数值矩阵。
    imuCfg.matVariable = '';
    imuCfg.rawTs = 0.002;
    imuCfg.samplesPerUpdate = round(ts/imuCfg.rawTs);
% 前三列为加速度计输出，后三列为陀螺输出，陀螺单位为 rad/s。
    imuCfg.accelCols = [1 2 3];
    imuCfg.gyroCols  = [4 5 6];
% 该 MAT 文件第二个加速度计通道静止均值接近 9.80，因此这里按 m/s^2 处理。
% 只有当文件中静止状态下的加速度数值约为 +/-1 时，才应改为 'g'。
    imuCfg.accelUnit = 'mps2';
% 从记录的物理输出中扣除的标定零偏。
    imuCfg.accelBiasMps2 = [0; 0; 0];
    imuCfg.gyroBiasRadPerSec = [0; 0; 0];
    imuCfg.autoEstimateGyroBias = true;
    imuCfg.gyroBiasStaticSamples = 10000;
% 数据坐标轴到导航平台坐标轴的转换关系：[x_data,y_data,z_data] -> [x_platform,z_platform,y_platform]。
% 这样可将数据第 2 轴上测得的正向重力映射到平台第 3 轴竖直方向。
    imuCfg.CdataToPlatform = [1 0 0; 0 0 -1; 0 1 0];
% 四环矩形路线仿真总时长按实际走完一圈设置。
% 一圈长度约 61.9 km，速度 15 m/s 时总时长约 1.15 h。
    route = make_fourth_ring_route(39.907456*arcdeg, 116.273987*arcdeg, 50, 15, 300);
    simTime = route.lapTime;
    route.simTime = simTime;
    route.simHours = simTime/3600;
    route.simLapCount = 1;
% 参考平台姿态四元数 qnpRef：单位四元数表示平台坐标系与导航系完全重合。
    [posRef, vnRef, attRef, trueEnu0] = fourth_ring_route_state(route, 0);
    qnpRef = [1; 0; 0; 0];
% qnp 是实际用于积分更新的平台姿态四元数，后面会叠加初始失准角。
    qnp = qnpRef;
% 把载体相对导航系的姿态角 attb0 转成四元数 qnb0。
    qnb0 = a2qua(attRef);
% Cbn0 是导航系速度投影到载体系/里程计坐标系的方向余弦矩阵，用作里程计观测矩阵。
% 车辆航向只用于里程计坐标分解，不改变平台姿态 qnp 的参考基准。
    Cbn0 = q2mat(qnb0)';
% 里程计相对于载体/车体系的安装误差角。
% 若设为 zeros(3,1)，则表示不考虑该安装误差。
    odoInstallErr = zeros(3,1);
    Cbo = rv2m(odoInstallErr);
    CbnOdo = Cbo*q2mat(a2qua(attRef))';
% IMU 相对于平台/载体系的安装误差角。
% 若设为 zeros(3,1)，则表示不考虑 IMU 安装误差。
    imuInstallErr = zeros(3,1);
    Cimu = rv2m(imuInstallErr);
% 初始参考位置：[纬度; 经度; 高度]，纬经度用弧度，高度用米。这里约为西安附近 34°N、108°E、高度 100 m。
% pos 是惯导解算位置，初值与参考位置相同，后续会因 IMU 误差产生漂移。
    pos = posRef;
% vn 是惯导解算速度，初值设为参考速度。
    vn = vnRef;
% 初始平台失准角 phi，单位为弧度：东向、北向为 0.1/0.2 角分，天向为 3 角分。
    phi = zeros(3,1);
% 把初始失准角注入平台姿态四元数，使后续误差曲线能体现初始姿态误差的收敛过程。
    qnp = qaddphi(qnp, phi);
    qnp0 = qnp;

% （2）滤波器参数初始化
% 软件接口建议：滤波器输入为当前惯导输出(qnp/vn/pos)、IMU周期平均比力、
% 以及外部速度/里程计观测；滤波器输出为 15 维误差状态 Xk。
% 里程计观测噪声标准差，分别对应载体系三个方向速度测量误差，单位 m/s。
    rk = [0.05; 0.02; 0.02];
% 观测噪声协方差 Rk，由标准差平方组成。
    Rk = diag(rk)^2;
% 初始状态协方差 P0：反映滤波器对初始误差的不确定性假设。
% 顺序依次为姿态误差、速度误差、位置误差、陀螺零偏、加速度计零偏。
    P0 = diag([[0.1; 0.1; 10]*arcdeg; [1; 1; 1]; [[10; 10]/Re; 10]; ...
               [0.1; 0.1; 0.1]*dph; [80; 90; 100]*ug])^2;
% 里程计观测模型 Hk：观测量是载体系速度残差 vb-odo。
% 对 15 状态而言，它主要直接观测速度误差项，所以中间 3 列是 Cbn0。
    Hk = [zeros(3), CbnOdo, zeros(3,9)];

% （3）惯性器件误差参数初始化
% 注意：这部分只用于“仿真 IMU 数据加误差”。如果将来采用真实惯导/IMU输出，
% 主控嵌入式软件一般不需要在这里人为设置 eb/web/db/wdb 来生成传感器误差。
% 陀螺误差参数：
%   eb  是常值零偏，单位由 deg/h 转成 rad/s；
%   web 是角随机游走强度，单位由 deg/sqrt(h) 转成 rad/sqrt(s)。
    eb  = [0.01; 0.015; 0.02]*dph;       web = [0.001; 0.001; 0.001]*dpsh;
% 加速度计误差参数：
%   db  是常值零偏，单位 ug；
%   wdb 是加速度随机噪声强度。
    db  = [80; 90; 100]*ug;              wdb = [1; 1; 1]*ugpsHz;

% 滤波器初始化收尾：Qk 使用上面定义的惯性器件随机误差强度。
% 过程噪声协方差 Qk：状态顺序为 [姿态误差3, 速度误差3, 位置误差3, 陀螺零偏3, 加速度计零偏3]。
% 这里只给陀螺随机游走和加速度计随机噪声建模，后 9 个状态的过程噪声置零。
    Qk = diag([web; wdb; zeros(9,1)])^2*nts;
% 初始化卡尔曼滤波器结构体，包括 Xk、Pk、Qk、Rk、状态转移矩阵和观测矩阵。
    kf = kfinit(Qk, Rk, P0, zeros(15), Hk);
    imuReader = [];
    if imuCfg.useRealData
        if strcmpi(imuCfg.source, 'mat')
            imuReader = initRealImuMatReader(imuCfg);
        else
            imuReader = initRealImuReader(imuCfg);
        end
        if strcmpi(imuCfg.source, 'mat')
            simTime = floor(size(imuReader.data,1)/imuCfg.samplesPerUpdate)*nts;
            if imuCfg.autoEstimateGyroBias
                biasSamples = min(imuCfg.gyroBiasStaticSamples, size(imuReader.data,1));
                gyroRaw = double(imuReader.data(1:biasSamples, imuCfg.gyroCols));
                gyroPlatform = (imuCfg.CdataToPlatform*gyroRaw.').';
                eth0 = earth(pos, vn);
                gyroExpectedPlatform = q2mat(qnp)'*eth0.wnie;
                imuCfg.gyroBiasRadPerSec = mean(gyroPlatform,1).' - gyroExpectedPlatform;
                imuReader.cfg = imuCfg;
                fprintf('Estimated gyro residual bias [rad/s]: [%.9g %.9g %.9g]\n', imuCfg.gyroBiasRadPerSec);
            end
        end
    end
% 根据总时长和更新周期计算循环次数。
    len = fix(simTime/nts);
% err 每行保存一次真实误差：[姿态误差3, 速度误差3, 位置误差3, 时间]。
    err = zeros(len, 10);
    nav = zeros(len, 11);
% xkpk 每行保存滤波状态估计 Xk、协方差对角线 diag(Pk) 和时间，用于后续画估计误差/收敛情况。
    xkpk = zeros(len, 2*kf.n+1);
    trajTime = zeros(len, 1);
    trajTrueEnu = zeros(len, 3);
    trajNavEnu = zeros(len, 3);
    trajTruePos = zeros(len, 3);
    trajNavPos = zeros(len, 3);
    trajTrueVn = zeros(len, 3);
% kk 是保存数组的行号，t 是当前仿真时间。
    kk = 1; t = 0;
% 主循环：每次循环推进 nts 秒，顺序为 生成 IMU → 惯导积分 → 卡尔曼预测/量测更新 → 反馈校正 → 保存误差。
    for k = 1:len
        t = t + nts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 模块二：惯性数据获取与惯性导航解算模块
% 输入：模块一初始化后的地理参数、参考状态、器件误差参数
% 输出：惯导速度 vn、姿态 qnp、位置 pos，以及本周期地球参数 eth
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% （1）惯性数据生成/获取
% 当前代码生成仿真 IMU 数据 wm1/vm1。
% 如果将来采用真实 IMU 输出：下面 ethRef/wm/vm/imuadderr 这一段由驱动或数据接口替换，
% 主控嵌入式软件通常直接读取真实角增量 wm1 和速度增量 vm1，不再生成理想 IMU 或人为加误差。
% 用参考位置和参考速度计算地球参数，作为理想 IMU 增量的来源。
        posRefPrev = posRef;
        vnRefPrev = vnRef;
        [posRefNext, vnRefNext, attRef, trueEnuNext] = fourth_ring_route_state(route, t);
        CbnOdo = Cbo*q2mat(a2qua(attRef))';
        ethRef = earth(posRefPrev, vnRefPrev);
        if imuCfg.useRealData
            if strcmpi(imuCfg.source, 'mat')
                [wm1, vm1, imuEof, imuReader] = readRealImuMatBlock(imuReader, imuCfg.samplesPerUpdate, imuCfg.rawTs, g0);
            else
                [wm1, vm1, imuEof] = readRealImuBlock(imuReader, imuCfg.samplesPerUpdate, imuCfg.rawTs, g0, arcdeg);
            end
            if imuEof
                fprintf('真实 IMU 数据在 %.3f s 处结束。\n', t-nts);
                break;
            end
        else
% 理想陀螺输出 wm：平台静止/匀速时，陀螺主要感受到导航系相对惯性系角速度 wnin。
% 每个小采样周期的角增量约为 wnin*ts，这里复制 nn 行供圆锥补偿函数使用。
        wm = repmat((ethRef.wnin*ts)', nn, 1);
% 理想加速度计输出 vm：静基座时主要抵消重力和运输项，写成 -gcc*ts 的速度增量。
        accRef = (vnRefNext - vnRefPrev)/nts;
        vm = repmat(((accRef - ethRef.gcc)*ts)', nn, 1);
% 给理想 IMU 增量加入陀螺零偏、陀螺随机游走、加速度计零偏和加速度随机噪声。
% 在加入传感器零偏和随机噪声之前，先施加 IMU 安装误差。
        wm = (Cimu*wm')';
        vm = (Cimu*vm')';
        [wm1, vm1] = imuadderr(wm, vm, eb, web, db, wdb, ts);
        end
% 与传感器输出同步生成或获取里程计输出。
% 里程计量测仍按 1 Hz 频率产生。
        if mod(t,1) < nts
            odo = CbnOdo*vnRefNext + rk.*randn(3,1);  % 里程计输出
        end

% （2）惯性导航解算
% 输入为本周期 IMU 角增量 wm1 和速度增量 vm1。
% 输出为当前姿态 qnp、速度 vn、位置 pos；这是后续滤波模块的主要输入。
% 平台惯导机械编排更新：由带误差的 IMU 增量积分得到新的姿态 qnp、速度 vn 和位置 pos。
        [qnp, vn, pos, eth] = pnsupdate(qnp, vn, pos, wm1, vm1, ts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 模块三：Kalman滤波、15维状态输出与仿真误差统计模块
% 输入：模块二输出的 qnp/vn/pos/eth、IMU增量 wm1/vm1、外部速度观测 odo
% 输出：15维滤波状态 kf.Xk、校正后的 qnp/vn/pos、仿真误差 err
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% （1）对速度、姿态、位置进行滤波处理，输出 15 维误差状态。
% 15维状态定义：
%   Xk( 1: 3) = 姿态误差 [phiE; phiN; phiU]
%   Xk( 4: 6) = 速度误差 [dvE; dvN; dvU]
%   Xk( 7: 9) = 位置误差 [dLat; dLon; dH]
%   Xk(10:12) = 陀螺零偏估计 [ebx; eby; ebz]
%   Xk(13:15) = 加速度计零偏估计 [dbx; dby; dbz]
% 构造离散状态转移矩阵 Phi ≈ I + F*nts。
% kff15 返回连续时间 15 状态误差模型矩阵 F，fb 取本周期平均比力。
        kf.Phikk_1 = eye(15) + kff15(eth, q2mat(qnp), sum(vm1,1)'/nts)*nts;
% 卡尔曼时间更新：根据状态转移矩阵传播状态估计和协方差。
        kf = kfupdate(kf);
% 里程计量测按 1 Hz 更新；由于导航周期是 0.2 s，这个条件大约每 5 次循环触发一次。
        if mod(t,1) < nts
% 把当前惯导估计速度从导航系投影到车体系，得到里程计应观测到的速度 vb。
            vb = CbnOdo*vn;
% odo 已在模块二中生成或获取，本段只使用该量测。
% 构造带噪声的里程计观测 odo：真实参考速度先投影到车体系，再叠加测量噪声。
% 本量测仍然只直接约束速度误差，所以 Hk 与初始化时相同。
            kf.Hk = [zeros(3), CbnOdo, zeros(3,9)];
% 卡尔曼量测更新：用速度残差 vb-odo 修正 15 维误差状态。
            kf = kfupdate(kf, vb-odo, 'M');
        end

% 15维滤波结果在 kf.Xk 中给出。嵌入式软件可根据闭环/开环方案决定是否执行下面的反馈校正。
% 闭环反馈校正：把滤波器估计出的姿态误差从平台姿态 qnp 中扣除，然后将该误差状态清零。
% 这样做相当于“估计一点、修正一点”，避免误差状态无限累积。
        qnp = qdelphi(qnp, kf.Xk(1:3)); kf.Xk(1:3) = 0;
% 速度闭环校正：用估计的速度误差修正惯导速度，并清零对应状态。
        vn  = vn - kf.Xk(4:6);          kf.Xk(4:6) = 0;
% 位置闭环校正：用估计的位置误差修正惯导位置，并清零对应状态。
        pos = pos - kf.Xk(7:9);         kf.Xk(7:9) = 0;

% （2）仿真误差计算
% 注意：这部分需要真实参考值 posRef/vnRef/qnpRef，所以只用于仿真验证。
% 如果采用真实数据，主控嵌入式软件没有真实值可相减，本段通常不进入实装代码。
% 经纬高误差换算成米：纬度误差乘地球半径，经度误差还要乘 cos(纬度)，高度误差直接相减。
        posRef = posRefNext;
        vnRef = vnRefNext;
        dpos = [(pos(1)-posRef(1))*Re;
                (pos(2)-posRef(2))*Re*cos(posRef(1));
                 pos(3)-posRef(3)];
% 保存真实误差：姿态误差由两个四元数求小角度误差，速度/位置直接与参考值相减。
        err(kk,:) = [qq2phi(qnp, qnpRef); vn-vnRef; dpos; t]';
        nav(kk,:) = [qnp; vn; pos; t]';
% 保存滤波器内部量：前 15 列为状态估计，后 15 列为协方差对角线，最后一列为时间。
        xkpk(kk,:) = [kf.Xk; diag(kf.Pk); t]';
        trueEnu = trueEnuNext;
        navEnu = trueEnu + [dpos(2); dpos(1); dpos(3)];
        trajTime(kk) = t;
        trajTrueEnu(kk,:) = trueEnu';
        trajNavEnu(kk,:) = navEnu';
        trajTruePos(kk,:) = posRef';
        trajNavPos(kk,:) = pos';
        trajTrueVn(kk,:) = vnRef';
% 参考位置按参考速度推进；静基座 vnRef=0 时参考位置不变，但保留该式可兼容匀速场景。
        kk = kk + 1;
% 每经过约 1 小时在命令行打印一次进度。当前 simTime=3600，所以结束时打印 1 h。
        if mod(t,3600) < nts
            fprintf('%s: %.0f h\n', sceneName, t/3600);
        end
    end
    if imuCfg.useRealData && strcmpi(imuCfg.source, 'text') && ~isempty(imuReader) && imuReader.fid > 0
        fclose(imuReader.fid);
    end
% 如果由于 kk 计数预分配多出尾部空行，这里裁掉未使用的部分。
    err(kk:end,:) = [];
    xkpk(kk:end,:) = [];
    nav(kk:end,:) = [];
    trajTime(kk:end,:) = [];
    trajTrueEnu(kk:end,:) = [];
    trajNavEnu(kk:end,:) = [];
    trajTruePos(kk:end,:) = [];
    trajNavPos(kk:end,:) = [];
    trajTrueVn(kk:end,:) = [];
% 将本次仿真的关键信息统一打包到 result 结构体，便于保存、画图和后续分析。
    result.sceneName = sceneName;
    result.err = err;
    result.nav = nav;
    result.nav_name = {'q0','q1','q2','q3','vE','vN','vU','lat','lon','h','t'};
    result.xkpk = xkpk;
    result.x15 = xkpk(:,1:15);
    result.x15_name = {'phiE','phiN','phiU','dvE','dvN','dvU', ...
                       'dLat','dLon','dH','ebx','eby','ebz','dbx','dby','dbz'};
    result.eb = eb;
    result.db = db;
    result.rk = rk;
    result.vnRef = vnRef;
    result.attb0 = attb0;
    result.qnp0 = qnp0;
    result.odoInstallErr = odoInstallErr;
    result.imuInstallErr = imuInstallErr;
    result.imuCfg = imuCfg;
    result.route = route;
    result.traj.time = trajTime;
    result.traj.trueEnu = trajTrueEnu;
    result.traj.navEnu = trajNavEnu;
    result.traj.truePos = trajTruePos;
    result.traj.navPos = trajNavPos;
    result.traj.trueVn = trajTrueVn;
% 计算各类误差的均方值，便于用一个数字评价姿态、速度、位置和零偏估计效果。
    if imuCfg.useRealData
        result.mse = [];
        fprintf('真实 IMU 滤波完成：%.3f s，最终速度 [%.6f %.6f %.6f] m/s。\n', ...
                nav(end,11), nav(end,5), nav(end,6), nav(end,7));
        plot_real_imu_navigation(result);
    else
        result.mse = calc_platform_odo_mse(err, xkpk, eb, db);
% 在命令行打印均方误差。
    print_platform_odo_mse(result);
% 绘制姿态误差、速度误差、位置误差、陀螺零偏估计和加速度计零偏估计曲线。
    plot_platform_odo_result(result);
    end
end

% ============================================================
% 误差统计函数：把时序误差转换成各通道均方误差 MSE
% ============================================================
function route = make_fourth_ring_route(lat0, lon0, h0, speed, turnRadius)
    width = 19000;
    height = 11000;
    startFromNorth = 0.70;
    startY = height*(1 - startFromNorth);
    arcLen = pi/2*turnRadius;
    segLen = [startY, arcLen, width, arcLen, height, arcLen, width, arcLen, height-startY];
    route.lat0 = lat0;
    route.lon0 = lon0;
    route.h0 = h0;
    route.pos0 = [lat0; lon0; h0];
    route.width = width;
    route.height = height;
    route.startY = startY;
    route.speed = speed;
    route.turnRadius = turnRadius;
    route.segLen = segLen;
    route.cumLen = [0, cumsum(segLen)];
    route.lapDistance = route.cumLen(end);
    route.lapTime = route.lapDistance/speed;
    route.wukesongLatDeg = lat0*180/pi;
    route.wukesongLonDeg = lon0*180/pi;
    route.disturb.lateralAmp = [0.8, 0.25, 0.08];
    route.disturb.lateralWave = [1200, 180, 37];
    route.disturb.lateralPhase = [0.4, 1.2, 2.0];
    route.disturb.alongAmp = [1.5, 0.3];
    route.disturb.alongWave = [2500, 260];
    route.disturb.alongPhase = [0.7, 0.1];
    route.disturb.diffStep = 0.2;
end

function [posRef, vnRef, attRef, enu] = fourth_ring_route_state(route, t)
    s = mod(route.speed*t, route.lapDistance);
    [enu, vnRef] = fourth_ring_actual_local_state(route, s);
    posRef = enu_to_pos(route.pos0, enu);
    yaw = atan2(vnRef(1), vnRef(2));
    attRef = [0; 0; yaw];
end

function [enu, vn] = fourth_ring_actual_local_state(route, s)
    enu = fourth_ring_actual_local_position(route, s);
    ds = route.disturb.diffStep;
    enuPlus = fourth_ring_actual_local_position(route, s + ds);
    enuMinus = fourth_ring_actual_local_position(route, s - ds);
    vn = route.speed*(enuPlus - enuMinus)/(2*ds);
end

function enu = fourth_ring_actual_local_position(route, s)
    s = mod(s, route.lapDistance);
    dAlong = harmonic_disturbance(s, route.disturb.alongAmp, ...
                                  route.disturb.alongWave, route.disturb.alongPhase);
    [baseEnu, baseVn] = fourth_ring_local_state(route, s + dAlong);
    speed2d = hypot(baseVn(1), baseVn(2));
    tangent = baseVn(1:2)/speed2d;
    normal = [-tangent(2); tangent(1)];
    dLat = harmonic_disturbance(s, route.disturb.lateralAmp, ...
                                route.disturb.lateralWave, route.disturb.lateralPhase);
    enu = baseEnu + [normal*dLat; 0];
end

function y = harmonic_disturbance(s, amp, wave, phase)
    y = 0;
    for ii = 1:numel(amp)
        y = y + amp(ii)*sin(2*pi*s/wave(ii) + phase(ii));
    end
end

function [enu, vn] = fourth_ring_local_state(route, s)
    s = mod(s, route.lapDistance);
    W = route.width;
    H = route.height;
    R = route.turnRadius;
    yStart = route.startY;
    v = route.speed;
    c = route.cumLen;

    if s < c(2)
        d = s - c(1);
        x = 0; y = yStart - d;
        vxy = [0; -v];
    elseif s < c(3)
        u = s - c(2);
        th = pi + u/R;
        x = R + R*cos(th); y = R*sin(th);
        vxy = v*[-sin(th); cos(th)];
    elseif s < c(4)
        d = s - c(3);
        x = R + d; y = -R;
        vxy = [v; 0];
    elseif s < c(5)
        u = s - c(4);
        th = -pi/2 + u/R;
        x = R + W + R*cos(th); y = R*sin(th);
        vxy = v*[-sin(th); cos(th)];
    elseif s < c(6)
        d = s - c(5);
        x = W + 2*R; y = d;
        vxy = [0; v];
    elseif s < c(7)
        u = s - c(6);
        th = u/R;
        x = W + R + R*cos(th); y = H + R*sin(th);
        vxy = v*[-sin(th); cos(th)];
    elseif s < c(8)
        d = s - c(7);
        x = W + R - d; y = H + R;
        vxy = [-v; 0];
    elseif s < c(9)
        u = s - c(8);
        th = pi/2 + u/R;
        x = R + R*cos(th); y = H + R*sin(th);
        vxy = v*[-sin(th); cos(th)];
    else
        d = s - c(9);
        x = 0; y = H - d;
        vxy = [0; -v];
    end

    enu = [x; y - yStart; 0];
    vn = [vxy(1); vxy(2); 0];
end

function pos = enu_to_pos(pos0, enu)
    eth0 = earth(pos0, [0; 0; 0]);
    lat = pos0(1) + enu(2)/eth0.RMh;
    lon = pos0(2) + enu(1)/eth0.clRNh;
    h = pos0(3) + enu(3);
    pos = [lat; lon; h];
end
function mse = calc_platform_odo_mse(err, xkpk, eb, db)
% arcmin/dph/ug 用于把弧度、rad/s、m/s^2 换算成更直观的角分、deg/h、ug。
    global arcmin dph ug
% 姿态误差前三列先由弧度换成角分，再分别对 E/N/U 三个通道求均方。
    mse.att_arcmin2 = mean((err(:,1:3)/arcmin).^2, 1);
% 速度误差本身单位就是 m/s，直接求均方。
    mse.vel_mps2 = mean(err(:,4:6).^2, 1);
% 位置误差已在主循环中换算为米，直接求均方。
    mse.pos_m2 = mean(err(:,7:9).^2, 1);
% 陀螺零偏估计误差：滤波估计值减真实零偏，再换算成 deg/h 后求均方。
    mse.gyro_bias_dph2 = mean(((xkpk(:,10:12)-eb')/dph).^2, 1);
% 加速度计零偏估计误差：滤波估计值减真实零偏，再换算成 ug 后求均方。
    mse.accel_bias_ug2 = mean(((xkpk(:,13:15)-db')/ug).^2, 1);
end

% ============================================================
% 打印统计结果：把 result.mse 中的各项指标按物理意义输出
% ============================================================
function print_platform_odo_mse(result)
    mse = result.mse;
% 每一行依次输出东、北、天或 x、y、z 三个通道的均方误差。
    fprintf('\n===== %s 四环一圈仿真的均方误差 =====\n', result.sceneName);
    fprintf('姿态误差 MSE [E N U] ((arcmin)^2):             %.6g  %.6g  %.6g\n', mse.att_arcmin2);
    fprintf('速度误差 MSE [E N U] ((m/s)^2):                %.6g  %.6g  %.6g\n', mse.vel_mps2);
    fprintf('位置误差 MSE [E N U] (m^2):                    %.6g  %.6g  %.6g\n', mse.pos_m2);
    fprintf('陀螺零偏误差 MSE [x y z] ((deg/h)^2):          %.6g  %.6g  %.6g\n', mse.gyro_bias_dph2);
    fprintf('加速度计零偏误差 MSE [x y z] ((ug)^2):         %.6g  %.6g  %.6g\n\n', mse.accel_bias_ug2);
end

% ============================================================
% 绘图函数：展示组合导航误差和零偏估计曲线
% ============================================================
function plot_platform_odo_result(result)
    global arcmin dph ug
% err 是真实误差时序；xkpk 是滤波器状态估计和协方差记录。
    err = result.err;
    xkpk = result.xkpk;
% err 最后一列为仿真时间，单位秒。
    tt = err(:,end);
% 如果仿真点数很多，最多抽取约 5000 个点绘图，避免图像刷新过慢。
    step = max(1, floor(numel(tt)/5000));
    idx = 1:step:numel(tt);
% 横轴换算为小时，更适合长时间导航误差展示。
    th = tt(idx)/3600;
% 子图 1：水平两个通道的平台姿态误差，单位角分。
    msplot(321, th, err(idx,1:2)/arcmin, 't / h', 'platform phi / arcmin');
    legend('phi E','phi N');
% 子图 2：天向/航向平台姿态误差，通常比水平姿态误差更难由里程计约束。
    msplot(322, th, err(idx,3)/arcmin, 't / h', 'platform phi U / arcmin');
% 子图 3：东、北、天三个速度误差。
    msplot(323, th, err(idx,4:6), 't / h', 'dvn / (m/s)');
    legend('dvE','dvN','dvU');
% 子图 4：东、北、天三个位置误差。
    msplot(324, th, err(idx,7:9), 't / h', 'dp / m');
    legend('dE','dN','dU');
% 子图 5：滤波器估计的三轴陀螺零偏，单位 deg/h。
    msplot(325, th, xkpk(idx,10:12)/dph, 't / h', 'gyro bias estimate / (deg/h)');
    legend('epsilon x','epsilon y','epsilon z');
% 子图 6：滤波器估计的三轴加速度计零偏，单位 ug。
    msplot(326, th, xkpk(idx,13:15)/ug, 't / h', 'accel bias estimate / ug');
    legend('nabla x','nabla y','nabla z');
    sgtitle([result.sceneName, ' platform INS/odometer one-day simulation']);
    plot_fourth_ring_route(result);
end

function plot_real_imu_navigation(result)
    global Re
    nav = result.nav;
    if isempty(nav)
        return;
    end

    step = max(1, floor(size(nav,1)/5000));
    idx = 1:step:size(nav,1);
    data = nav(idx,:);
    timeMin = data(:,11)/60;

    attDeg = zeros(numel(idx),3);
    for i = 1:numel(idx)
        qRelative = qmul(qconj(result.qnp0), data(i,1:4).');
        attDeg(i,:) = q2att321(qRelative)'*180/pi;
    end

    lat0 = nav(1,8);
    lon0 = nav(1,9);
    h0 = nav(1,10);
    posENU = [(data(:,9)-lon0)*Re*cos(lat0), ...
              (data(:,8)-lat0)*Re, ...
               data(:,10)-h0];

    figure('Name', 'Real IMU INS/Odometer Filtering', 'NumberTitle', 'off');
    subplot(4,1,1);
    plot(timeMin, attDeg, 'LineWidth', 1); grid on;
    xlabel('t / min'); ylabel('attitude change / deg');
    legend('roll change','pitch change','yaw change', 'Location', 'best');

    subplot(4,1,2);
    plot(timeMin, data(:,5:7), 'LineWidth', 1); grid on;
    xlabel('t / min'); ylabel('velocity / (m/s)');
    legend('vE','vN','vU', 'Location', 'best');

    subplot(4,1,3);
    plot(timeMin, posENU, 'LineWidth', 1); grid on;
    xlabel('t / min'); ylabel('relative position / m');
    legend('E','N','U', 'Location', 'best');

    subplot(4,1,4);
    plot(posENU(:,1), posENU(:,2), 'LineWidth', 1); grid on; axis equal;
    xlabel('East / m'); ylabel('North / m');
    title('Horizontal trajectory relative to the initial position');
    sgtitle([result.sceneName, ' real IMU INS/odometer filtering']);
end

function plot_fourth_ring_route(result)
    route = result.route;
    trueEnu = result.traj.trueEnu;
    navEnu = result.traj.navEnu;
    if isempty(trueEnu)
        return;
    end
    step = max(1, floor(size(trueEnu,1)/3000));
    idx = 1:step:size(trueEnu,1);

    figure('Name', 'Actual vehicle fourth-ring trajectory', 'NumberTitle', 'off');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile;
    plot(trueEnu(idx,1)/1000, trueEnu(idx,2)/1000, 'b-', 'LineWidth', 2.0); hold on;
    plot(navEnu(idx,1)/1000, navEnu(idx,2)/1000, 'r--', 'LineWidth', 1.0);
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    quiver(0, 0, 0, -0.8, 0, 'k', 'LineWidth', 1.2, 'MaxHeadSize', 0.8);
    grid on; axis equal;
    xlabel('East / km');
    ylabel('North / km');
    title('Full actual vehicle trajectory');
    legend('Actual vehicle trajectory', 'INS/odometer navigation trajectory', ...
           'Start / Wukesong anchor', 'Initial direction', 'Location', 'best');
    text(0.2, 0.2, sprintf('Wukesong %.6f N, %.6f E', ...
         route.wukesongLatDeg, route.wukesongLonDeg), 'FontSize', 9);

    zoomDistance = 1500;
    zoomIdx = find(result.traj.time*route.speed <= zoomDistance);
    nominalEnu = zeros(numel(zoomIdx), 3);
    for ii = 1:numel(zoomIdx)
        s = mod(route.speed*result.traj.time(zoomIdx(ii)), route.lapDistance);
        nominalEnu(ii,:) = fourth_ring_local_state(route, s)';
    end
    zoomTrue = trueEnu(zoomIdx,:) - nominalEnu(1,:);
    zoomNominal = nominalEnu - nominalEnu(1,:);

    nexttile;
    plot(zoomNominal(:,1), zoomNominal(:,2), 'k:', 'LineWidth', 1.2); hold on;
    plot(zoomTrue(:,1), zoomTrue(:,2), 'b-', 'LineWidth', 1.6);
    grid on; xlim([-2, 2]);
    xlabel('East / m');
    ylabel('North / m');
    title('Zoomed actual motion disturbance');
    legend('Ideal centerline', 'Actual disturbed trajectory', 'Location', 'best');
    totalHours = result.traj.time(end)/3600;
    lapCount = result.traj.time(end)*route.speed/route.lapDistance;
    sgtitle(sprintf('Fourth-ring vehicle trajectory, R=%.0f m, v=%.1f m/s, %.1f h / %.2f laps', ...
            route.turnRadius, route.speed, totalHours, lapCount));
    exportgraphics(gcf, 'real_imu_interface_fourth_ring_actual_vehicle_trajectory.png', 'Resolution', 200);
end
function att = q2att321(q)
    Cnb = q2mat(q);
    roll = atan2(Cnb(3,2), Cnb(3,3));
    pitch = asin(max(-1, min(1, -Cnb(3,1))));
    yaw = atan2(Cnb(2,1), Cnb(1,1));
    att = [roll; pitch; yaw];
end

% ============================================================
% 姿态角转四元数：att = [roll; pitch; yaw]，单位弧度
% 四元数排列采用 [q0; q1; q2; q3]，其中 q0 为标量部。
% ============================================================
function qnb = a2qua(att)
% 强制转换成列向量，避免输入为行向量时后续索引或矩阵运算出错。
    att = att(:);
% 欧拉角转四元数公式中需要用半角的 sin/cos。
    s = sin(att/2); c = cos(att/2);
    si = s(1); sj = s(2); sk = s(3);
    ci = c(1); cj = c(2); ck = c(3);
% 按 3-2-1 或本程序约定的姿态角顺序组合四元数分量。
    qnb = [ ci*cj*ck + si*sj*sk;
            si*cj*ck - ci*sj*sk;
            ci*sj*ck + si*cj*sk;
            ci*cj*sk + si*sj*ck ];
% 归一化可以抑制数值舍入导致的四元数模长偏离 1。
    qnb = qnormlz(qnb);
end

% ============================================================
% 反对称矩阵函数：把向量 v 构造成叉乘矩阵 [v×]
% 满足 askew(v)*a = cross(v,a)，常用于姿态误差和旋转向量公式。
% ============================================================
function m = askew(v)
% 保证 v 是 3×1 列向量。
    v = v(:);
% 该矩阵是反对称矩阵，对角线为 0，上下三角互为相反数。
    m = [ 0,    -v(3),  v(2);
          v(3),  0,    -v(1);
         -v(2),  v(1),  0    ];
end

% ============================================================
% 圆锥与划桨补偿 cnscul
% 输入 wm/vm 是一个导航更新周期内的若干个陀螺角增量和加速度计速度增量。
% 输出 phim 为合成旋转矢量，dvbm 为补偿后的速度增量。
% ============================================================
function [phim, dvbm] = cnscul(wm, vm)
% cs 是多子样圆锥/划桨补偿系数表；n 不同时选用对应行。
    cs = [ [2,    0,    0,    0,     0]/3;
           [9,   27,    0,    0,     0]/20;
           [54,  92,  214,    0,     0]/105;
           [250, 525, 650, 1375,     0]/504;
           [2315,4558,7296,7834,15797]/4620 ];
% wmm/vmm 是本周期所有小采样角增量、速度增量的直接累加。
    wmm = sum(wm,1);
    vmm = sum(vm,1);
% dphim 是圆锥补偿项，scullm 是划桨补偿项，先初始化为零。
    dphim = zeros(1,3);
    scullm = zeros(1,3);
% n 为本周期内的小采样个数。若 n=1，没有多子样补偿项。
    n = size(wm,1);
    if n > 1
% 对前 n-1 个子样加权，形成补偿所需的历史角增量/速度增量。
        csw = cs(n-1,1:n-1)*wm(1:n-1,:);
        csv = cs(n-1,1:n-1)*vm(1:n-1,:);
% 圆锥补偿：考虑不同子样角增量不共线时的二阶旋转误差。
        dphim  = cross(csw, wm(n,:));
% 划桨补偿：考虑旋转与速度增量耦合导致的二阶速度误差。
        scullm = cross(csw, vm(n,:)) + cross(csv, wm(n,:));
    end
% 合成后的姿态旋转矢量。
    phim = (wmm + dphim)';
% 合成后的速度增量，包含 0.5*cross(wmm,vmm) 的旋转速度耦合项和划桨补偿项。
    dvbm = (vmm + 0.5*cross(wmm,vmm) + scullm)';
end

% ============================================================
% 地球参数计算 earth
% 输入 pos=[纬度;经度;高度]、vn=[vE;vN;vU]，输出曲率半径、自转项、运输项和重力等。
% 这些量是惯导机械编排和误差方程的基础。
% ============================================================
function eth = earth(pos, vn)
% Re 为地球长半轴，ff 为扁率，wie 为地球自转角速度，g0 为赤道标准重力。
    global Re ff wie g0
    pos = pos(:); vn = vn(:);
% 由扁率计算第一偏心率，用于子午圈/卯酉圈曲率半径。
    ee = sqrt(2*ff-ff^2);
    e2 = ee^2;
% 保存纬度的 sin/cos/tan，后续多处重复使用。
    eth.sl = sin(pos(1));
    eth.cl = cos(pos(1));
    eth.tl = eth.sl/eth.cl;
    eth.sl2 = eth.sl*eth.sl;
    sl4 = eth.sl2*eth.sl2;
% sq = 1 - e^2 sin^2(L)，是椭球曲率半径公式中的公共项。
    sq = 1 - e2*eth.sl2;
    sq2 = sqrt(sq);
% RMh：子午圈曲率半径加高度，用于纬度更新和北向运输角速度。
    eth.RMh = Re*(1-e2)/sq/sq2 + pos(3);
% RNh：卯酉圈曲率半径加高度，用于经度更新和东向运输角速度。
    eth.RNh = Re/sq2 + pos(3);
% clRNh = cos(L)*(RN+h)，经度变化率分母会用到。
    eth.clRNh = eth.cl*eth.RNh;
% wnie：地球自转角速度在导航系 n 下的投影。
    eth.wnie = wie*[0; eth.cl; eth.sl];
    eth.pos = pos;
    eth.vn = vn;
% wnen：导航系相对地球系的运输角速度，由载体在地球表面的运动产生。
    eth.wnen = [-vn(2)/eth.RMh;
                 vn(1)/eth.RNh;
                 vn(1)/eth.RNh*eth.tl];
% wnin：导航系相对惯性系角速度 = 地球自转项 + 运输项。
    eth.wnin = eth.wnie + eth.wnen;
% wnien 在速度方程中用于构造哥氏/运输相关项。
    eth.wnien = eth.wnie + eth.wnin;
% 正常重力模型：随纬度和高度变化，纬度越高重力通常越大，高度越高重力越小。
    gLh = g0*(1 + 5.27094e-3*eth.sl2 + 2.32718e-5*sl4) - 3.086e-6*pos(3);
% 导航系采用东-北-天，重力沿天向负方向。
    eth.gn = [0; 0; -gLh];
% gcc 是速度更新中使用的等效重力/有害加速度项，包含重力和旋转引起的修正。
    eth.gcc = eth.gn - cross(eth.wnien, vn);
end

% ============================================================
% IMU 误差注入函数 imuadderr
% 在理想角增量 wm 和速度增量 vm 中加入常值零偏与白噪声/随机游走。
% ============================================================
function [wm, vm] = imuadderr(wm, vm, eb, web, db, wdb, ts)
% m 是当前更新周期内的子样数。
    m = size(wm,1);
% 离散白噪声标准差与 sqrt(ts) 成比例，因此先计算 sqrt(ts)。
    sts = sqrt(ts);
% 陀螺角增量误差 = 零偏*ts + 随机游走强度*sqrt(ts)*N(0,1)。
    wm = wm + [ ts*eb(1) + sts*web(1)*randn(m,1), ...
                ts*eb(2) + sts*web(2)*randn(m,1), ...
                ts*eb(3) + sts*web(3)*randn(m,1) ];
% 加速度计速度增量误差 = 加速度零偏*ts + 噪声强度*sqrt(ts)*N(0,1)。
    vm = vm + [ ts*db(1) + sts*wdb(1)*randn(m,1), ...
                ts*db(2) + sts*wdb(2)*randn(m,1), ...
                ts*db(3) + sts*wdb(3)*randn(m,1) ];
end

% ============================================================
% 平台惯导机械编排更新 pnsupdate
% 根据 IMU 增量更新平台姿态四元数 qnp、速度 vn、位置 pos。
% qnp 表示平台坐标系 p 到导航系 n 的姿态四元数。
% ============================================================
function [qnp, vn, pos, eth] = pnsupdate(qnp, vn, pos, wm, vm, ts)
% 当前周期内包含的 IMU 子样数。
    nn = size(wm,1);
    nts = nn*ts;
% 对子样 IMU 增量做圆锥/划桨补偿，得到本周期合成旋转矢量和速度增量。
    [phim, dvpm] = cnscul(wm, vm);
% 根据当前估计位置和速度计算地球参数。
    eth = earth(pos, vn);
% 速度预更新：
%   qmulv(qnp,dvpm) 把平台系速度增量转到导航系；
%   rv2m(-wnin*nts/2) 做半周期旋转补偿；
%   eth.gcc*nts 加入重力和旋转修正。
    vn1 = vn + rv2m(-eth.wnin*nts/2)*qmulv(qnp, dvpm) + eth.gcc*nts;
% 用梯形思想取更新前后速度平均值，用于位置更新。
    vn = (vn + vn1)/2;
% 位置更新：纬度由北向速度决定，经度由东向速度决定，高度由天向速度决定。
    pos = pos + [vn(2)/eth.RMh; vn(1)/eth.clRNh; vn(3)]*nts;
% 位置更新后，把速度正式替换为本周期末速度。
    vn = vn1;
% 姿态更新：左乘导航系旋转补偿，右乘 IMU 测得的平台旋转增量。
    qnp = qmul(rv2q(-eth.wnin*nts), qmul(qnp, rv2q(phim)));
% 四元数积分后归一化，避免长时间仿真中数值漂移。
    qnp = qnormlz(qnp);
end

% ============================================================
% 15 状态惯导误差方程矩阵 kff15
% 状态顺序：phi(姿态误差3), dvn(速度误差3), dpos(位置误差3), eb(陀螺零偏3), db(加速度计零偏3)。
% 输出 Ft 是连续时间误差状态矩阵 F，主程序中用 I+F*nts 近似离散化。
% ============================================================
function Ft = kff15(eth, Cnb, fb)
    global g0 ff
% tl/secL 等纬度相关量，以及经纬高和曲率半径倒数，后续用于构造误差矩阵。
    tl = eth.tl; secl = 1/eth.cl;
    L = eth.pos(1); h = eth.pos(3);
    f_RMh = 1/eth.RMh; f_RNh = 1/eth.RNh; f_clRNh = 1/eth.clRNh;
    f_RMh2 = f_RMh*f_RMh; f_RNh2 = f_RNh*f_RNh;
    vE_clRNh = eth.vn(1)*f_clRNh; vE_RNh2 = eth.vn(1)*f_RNh2;
    vN_RMh2 = eth.vn(2)*f_RMh2;
% Mp1：位置误差对地球自转投影项的影响。
    Mp1 = [0,              0, 0;
          -eth.wnie(3),    0, 0;
           eth.wnie(2),    0, 0];
% Mp2：位置误差对运输角速度项的影响，与速度和曲率半径有关。
    Mp2 = [0,              0,  vN_RMh2;
           0,              0, -vE_RNh2;
           vE_clRNh*secl,  0, -vE_RNh2*tl];
% beta 系数用于重力随纬度和高度变化的线性化近似。
    beta = 5.27094e-3; beta1 = (2*beta+ff)*ff/8; beta2 = 3.086e-6; beta3 = 8.08e-9;
% Mp3：位置误差对重力项的影响，即重力模型线性化。
    Mp3 = [0, 0, 0;
          -2*beta3*h, 0, -beta3*sin(2*L);
          -g0*(beta-4*beta1*cos(2*L))*sin(2*L), 0, beta2];
% Maa：姿态误差自身传播项，由导航系相对惯性系角速度决定。
    Maa = askew(-eth.wnin);
% Mav：速度误差对姿态误差的耦合，来自运输角速度对速度的依赖。
    Mav = [0,        -f_RMh, 0;
           f_RNh,     0,    0;
           f_RNh*tl,  0,    0];
% Map：位置误差对姿态误差的耦合，由地球自转投影和运输项共同组成。
    Map = Mp1 + Mp2;
% Mva：姿态误差对速度误差的影响，本质是比力方向错误导致的速度误差。
    Mva = askew(Cnb*fb(:));
% Mvv：速度误差自身传播项，包含旋转项和运输角速度相关项。
    Mvv = askew(eth.vn)*Mav - askew(eth.wnien);
% Mvp：位置误差对速度误差的影响，包含运输项变化和重力变化。
    Mvp = askew(eth.vn)*(Mp1+Map) + Mp3;
% Mpv：速度误差对位置误差的影响，对应经纬高运动学方程的线性化。
    Mpv = [0,        f_RMh, 0;
           f_clRNh,  0,    0;
           0,        0,    1];
% Mpp：位置误差自身传播项，由曲率半径和纬度变化引起。
    Mpp = [0,             0, -vN_RMh2;
           vE_clRNh*tl,   0, -vE_RNh2*secl;
           0,             0,  0];
% O33 是 3×3 零矩阵，用于拼接分块矩阵。
    O33 = zeros(3);
% 拼接完整 15×15 连续误差矩阵：
%   第1行块：姿态误差方程，陀螺零偏通过 -Cnb 进入姿态误差；
%   第2行块：速度误差方程，加速度计零偏通过 Cnb 进入速度误差；
%   第3行块：位置误差方程；
%   最后6行：本程序把陀螺/加速度计零偏建模为常值，所以对应导数为 0。
    Ft = [Maa, Mav, Map, -Cnb, O33;
          Mva, Mvv, Mvp,  O33, Cnb;
          O33, Mpv, Mpp,  O33, O33;
          zeros(6,15)];
end

% ============================================================
% 卡尔曼滤波器初始化 kfinit
% 把滤波所需矩阵和状态统一放入结构体 kf，便于主循环调用。
% ============================================================
function kf = kfinit(Qk, Rk, P0, Phikk_1, Hk, Gammak)
% Hk 的行数是量测维数 m，列数是状态维数 n。当前 m=3，n=15。
    [kf.m, kf.n] = size(Hk);
% 保存过程噪声协方差 Qk。
    kf.Qk = Qk;
% 保存观测噪声协方差 Rk。
    kf.Rk = Rk;
% Pk 是当前状态估计误差协方差，初值为 P0。
    kf.Pk = P0;
% Xk 是当前状态估计，初始认为误差估计为 0。
    kf.Xk = zeros(kf.n,1);
% Phikk_1 是离散状态转移矩阵，主循环中每步更新。
    kf.Phikk_1 = Phikk_1;
% Hk 是量测矩阵，里程计更新时也可重新赋值。
    kf.Hk = Hk;
% 若没有显式给噪声驱动矩阵 Gamma，默认过程噪声直接作用在所有状态上。
    if nargin < 6
        kf.Gammak = eye(kf.n);
    else
        kf.Gammak = Gammak;
    end
end

% ============================================================
% 卡尔曼滤波更新 kfupdate
% TimeMeasBoth 控制更新类型：'T' 只时间更新，'M' 只量测更新，'B' 时间+量测同时更新。
% ============================================================
function kf = kfupdate(kf, Zk, TimeMeasBoth)
% 只传入 kf 时，默认执行时间更新。
    if nargin == 1
        TimeMeasBoth = 'T';
% 传入 kf 和 Zk 但不指定类型时，默认执行时间更新+量测更新。
    elseif nargin == 2
        TimeMeasBoth = 'B';
    end
% 时间更新：用状态转移矩阵传播状态估计和协方差。
    if TimeMeasBoth == 'T' || TimeMeasBoth == 'B'
% 状态预测：X(k|k-1)=Phi*X(k-1|k-1)。
        kf.Xkk_1 = kf.Phikk_1*kf.Xk;
% 协方差预测：P(k|k-1)=Phi*P*Phi' + Gamma*Q*Gamma'。
        kf.Pkk_1 = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';
    else
        kf.Xkk_1 = kf.Xk;
        kf.Pkk_1 = kf.Pk;
    end
% 量测更新：有新观测时计算卡尔曼增益并修正状态。
    if TimeMeasBoth == 'M' || TimeMeasBoth == 'B'
% 状态与量测的互协方差 P_xz = P*H'。
        kf.PXZkk_1 = kf.Pkk_1*kf.Hk';
% 量测预测协方差 P_zz = H*P*H' + R。
        kf.PZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;
% 卡尔曼增益 K = P_xz / P_zz，MATLAB 右除等价于乘 P_zz 的逆。
        kf.Kk = kf.PXZkk_1/kf.PZkk_1;
% 状态修正：预测值加上 K 乘以量测残差。
        kf.Xk = kf.Xkk_1 + kf.Kk*(Zk - kf.Hk*kf.Xkk_1);
% 协方差修正：这里采用简化形式 P = P - K*Pzz*K'。
        kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZkk_1*kf.Kk';
    else
        kf.Xk = kf.Xkk_1;
        kf.Pk = kf.Pkk_1;
    end
% 数值对称化，防止浮点误差让协方差矩阵出现轻微非对称。
    kf.Pk = (kf.Pk + kf.Pk')/2;
end

% ============================================================
% 绘图辅助函数 msplot
% mnp 采用 MATLAB subplot 的三位数写法，例如 321 表示 3行2列第1幅图。
% ============================================================
function msplot(mnp, x, y, xstr, ystr)
% 当子图编号个位为 1 时，新开一个 figure，便于一组子图放在同一窗口。
    if mod(mnp,10) == 1
        figure;
    end
% 绘制曲线并打开网格。
    subplot(mnp); plot(x, y); grid on;
% 如果只给了 xstr/ystr 中的一个标签，则默认横轴为 t/s。
    if nargin == 4
        ystr = xstr;
        xstr = 't / s';
    end
% 设置坐标轴标签。
    xlabel(xstr); ylabel(ystr);
end

% ============================================================
% 四元数转方向余弦矩阵 q2mat
% 输入 qnb 为导航相关四元数，输出 Cnb，用于向量坐标变换。
% ============================================================

function q = m2qua(C)
    tr = trace(C);
    if tr > 0
        s = 2*sqrt(tr + 1);
        q = [0.25*s; (C(3,2)-C(2,3))/s; (C(1,3)-C(3,1))/s; (C(2,1)-C(1,2))/s];
    elseif C(1,1) > C(2,2) && C(1,1) > C(3,3)
        s = 2*sqrt(1 + C(1,1) - C(2,2) - C(3,3));
        q = [(C(3,2)-C(2,3))/s; 0.25*s; (C(1,2)+C(2,1))/s; (C(1,3)+C(3,1))/s];
    elseif C(2,2) > C(3,3)
        s = 2*sqrt(1 + C(2,2) - C(1,1) - C(3,3));
        q = [(C(1,3)-C(3,1))/s; (C(1,2)+C(2,1))/s; 0.25*s; (C(2,3)+C(3,2))/s];
    else
        s = 2*sqrt(1 + C(3,3) - C(1,1) - C(2,2));
        q = [(C(2,1)-C(1,2))/s; (C(1,3)+C(3,1))/s; (C(2,3)+C(3,2))/s; 0.25*s];
    end
    q = qnormlz(q);
end
function Cnb = q2mat(qnb)
% 使用前先归一化，保证方向余弦矩阵近似正交。
    qnb = qnormlz(qnb(:));
% 预先计算四元数分量乘积，减少重复运算并使矩阵公式更清楚。
    q11 = qnb(1)*qnb(1); q12 = qnb(1)*qnb(2); q13 = qnb(1)*qnb(3); q14 = qnb(1)*qnb(4);
    q22 = qnb(2)*qnb(2); q23 = qnb(2)*qnb(3); q24 = qnb(2)*qnb(4);
    q33 = qnb(3)*qnb(3); q34 = qnb(3)*qnb(4);
    q44 = qnb(4)*qnb(4);
% 按四元数到方向余弦矩阵的标准公式构造 3×3 矩阵。
    Cnb = [ q11+q22-q33-q44,  2*(q23-q14),      2*(q24+q13);
            2*(q23+q14),      q11-q22+q33-q44,  2*(q34-q12);
            2*(q24-q13),      2*(q34+q12),      q11-q22-q33+q44 ];
end

% ============================================================
% qaddphi：给姿态四元数加入小失准角 phi
% 在主程序中用于模拟平台初始姿态误差。
% ============================================================
function qpb = qaddphi(qnb, phi)
% 失准角先转成旋转四元数，再左乘到原姿态四元数上。
    qpb = qmul(rv2q(-phi(:)), qnb(:));
    qpb = qnormlz(qpb);
end

% ============================================================
% 四元数共轭 qconj：标量部不变，向量部取反
% 对单位四元数，共轭等于逆，用于反向旋转或构造误差四元数。
% ============================================================
function qout = qconj(qin)
    qin = qin(:);
    qout = [qin(1); -qin(2:4)];
end

% ============================================================
% qdelphi：从姿态四元数中扣除小失准角 phi
% 在闭环滤波反馈中用于姿态校正。
% ============================================================
function qnb = qdelphi(qpb, phi)
% 与 qaddphi 符号相反，用估计的小角度误差修正当前平台姿态。
    qnb = qmul(rv2q(phi(:)), qpb(:));
    qnb = qnormlz(qnb);
end

% ============================================================
% 四元数乘法 qmul：q = q1*q2
% 姿态连续旋转时用四元数乘法进行组合。
% ============================================================
function q = qmul(q1, q2)
% 保证输入均为列向量，避免行列维度不一致。
    q1 = q1(:); q2 = q2(:);
% 第一行是标量部，后三行是向量部，按 Hamilton 四元数乘法展开。
    q = [ q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);
          q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3);
          q1(1)*q2(3) + q1(3)*q2(1) + q1(4)*q2(2) - q1(2)*q2(4);
          q1(1)*q2(4) + q1(4)*q2(1) + q1(2)*q2(3) - q1(3)*q2(2) ];
end

% ============================================================
% qmulv：用四元数旋转三维向量
% 计算形式为 q * [0;v] * q^*，输出结果取向量部。
% ============================================================
function vo = qmulv(q, vi)
    vi = vi(:);
% 把三维向量嵌入纯四元数，标量部为 0。
    qi = [0; vi];
% 左乘 q、右乘 q 的共轭，完成坐标旋转。
    qo = qmul(qmul(q(:), qi), qconj(q(:)));
    vo = qo(2:4);
end

% ============================================================
% 四元数归一化 qnormlz
% 长时间积分时必须不断归一化，否则四元数模长会因数值误差偏离 1。
% ============================================================
function qnb = qnormlz(qnb)
    qnb = qnb(:);
% nm 是四元数模长平方。
    nm = qnb'*qnb;
% 若模长过小，说明数值异常，直接重置为单位四元数以避免除零。
    if nm < 1e-6
        qnb = [1; 0; 0; 0];
    else
        qnb = qnb/sqrt(nm);
    end
end

% ============================================================
% qq2phi：由两个四元数计算小姿态误差角
% qpb 是当前平台姿态，qnb 是参考姿态，输出 phi 为旋转矢量形式的小角度误差。
% ============================================================
function phi = qq2phi(qpb, qnb)
% 先构造误差四元数 qerr，再转成旋转矢量。
    qerr = qmul(qnb(:), qconj(qpb(:)));
    phi = q2rv(qerr);
end

% ============================================================
% 四元数转旋转矢量 q2rv
% 对小角度时，旋转矢量近似等于 2*q 的向量部。
% ============================================================
function rv = q2rv(q)
% 归一化并统一四元数符号，避免同一姿态由 q 和 -q 两种形式导致跳变。
    q = qnormlz(q(:));
    if q(1) < 0
        q = -q;
    end
% nmhalf 是旋转角的一半，通过 acos(q0) 得到，并用 max/min 防止数值越界。
    nmhalf = acos(max(min(q(1),1),-1));
% 大角度用精确比例系数，小角度用极限值 2，避免 sin(nmhalf) 太小导致数值不稳定。
    if nmhalf > 1e-20
        b = 2*nmhalf/sin(nmhalf);
    else
        b = 2;
    end
    rv = b*q(2:4);
end

% ============================================================
% 旋转矢量转方向余弦矩阵 rv2m
% 使用 Rodrigues 公式：C = I + a*[rv×] + b*[rv×]^2。
% ============================================================
function m = rv2m(rv)
    rv = rv(:);
% nm2 是旋转矢量模长平方。
    nm2 = rv'*rv;
% 小角度时用泰勒展开计算 a、b，避免 sin(x)/x 和 (1-cos(x))/x^2 的数值问题。
    if nm2 < 1e-8
        a = 1 - nm2*(1/6 - nm2/120);
        b = 0.5 - nm2*(1/24 - nm2/720);
    else
        nm = sqrt(nm2);
        a = sin(nm)/nm;
        b = (1-cos(nm))/nm2;
    end
% VX 是旋转矢量对应的叉乘矩阵。
    VX = askew(rv);
% Rodrigues 公式得到方向余弦矩阵。
    m = eye(3) + a*VX + b*VX^2;
end

% ============================================================
% 旋转矢量转四元数 rv2q
% 用于把小角度误差或 IMU 合成旋转矢量转换为姿态更新四元数。
% ============================================================
function q = rv2q(rv)
    rv = rv(:);
% 旋转矢量模长平方。
    nm2 = rv'*rv;
% 小角度时对 cos(|rv|/2) 和 sin(|rv|/2)/|rv| 使用泰勒展开。
    if nm2 < 1.0e-8
        q0 = 1 - nm2*(1/8 - nm2/384);
        s  = 1/2 - nm2*(1/48 - nm2/3840);
    else
% 一般角度下使用精确三角函数计算。
        nm = sqrt(nm2);
        q0 = cos(nm/2);
        s  = sin(nm/2)/nm;
    end
% 拼成 [标量部; 向量部] 后归一化。
    q = qnormlz([q0; s*rv]);
end

% ============================================================
% 全局常量函数 gvar
% 把惯导仿真中常用的地球参数、单位换算和绘图线型放入 global 变量。
% ============================================================
function gvar
% 这些 global 变量被 earth、主仿真、误差统计和绘图等多个局部函数共享。
global GM Re ff wie ge gp g0 ug arcdeg arcmin arcsec hur dph dpsh ugpsHz lsc
% 地球引力常数 GM，当前代码中保留但未直接使用。
GM  = 3.986004415e14;
% Re：WGS-84 地球长半轴/赤道半径，单位 m。
Re  = 6.378136998405e6;
% wie：地球自转角速度，单位 rad/s。
wie = 7.2921151467e-5;
% ff：WGS-84 扁率。
ff  = 1/298.257223563;
% ee/e2/Rp 在本函数中计算后未被其他函数直接使用，#ok<NASGU> 用于抑制 MATLAB 未使用变量警告。
ee  = sqrt(2*ff-ff^2); 
e2  = ee^2;            
Rp  = (1-ff)*Re;       
% gE/ge/gp/g0：赤道、极区和本程序采用的标准重力参数。
gE = 9.780325333434361; %#ok<NASGU>
ge = gE;
gp = 9.832184935381024;
g0 = ge;
% ug：微 g 单位，1 ug = g0*1e-6，用于加速度计零偏换算。
ug = g0*1e-6;
% 角度单位换算：度、角分、角秒到弧度。
arcdeg = pi/180;
arcmin = arcdeg/60;
arcsec = arcmin/60;
% hur 表示 1 小时的秒数。
hur    = 3600;
% dph：deg/hour 换算到 rad/s，用于陀螺常值零偏。
dph    = arcdeg/hur;
% dpsh：deg/sqrt(hour) 换算到 rad/sqrt(second)，用于陀螺随机游走。
dpsh   = arcdeg/sqrt(hur);
% ugpsHz：加速度噪声单位，这里相当于 ug/sqrt(Hz) 的尺度。
ugpsHz = ug/sqrt(1);
% lsc：绘图线型和颜色列表，当前 msplot 中没有直接使用，但保留以兼容其他脚本风格。
lsc = {'-k','-b','-r','-m','-g','--k','--b','--r','--m','--g', ...
       ':k',':b',':r',':m',':g'};
end

% ============================================================
% 真实 IMU 文本文件接口
% ============================================================
function reader = initRealImuMatReader(cfg)
    fileInfo = whos('-file', cfg.file);
    if isempty(fileInfo)
        error('真实 IMU 的 MAT 文件中没有找到任何变量：%s', cfg.file);
    end

    variableName = cfg.matVariable;
    if isempty(variableName)
        neededCols = max([cfg.accelCols(:); cfg.gyroCols(:)]);
        variableName = '';
        for i = 1:numel(fileInfo)
            sizeInfo = fileInfo(i).size;
            if numel(sizeInfo) == 2 && (sizeInfo(2) >= neededCols || sizeInfo(1) >= neededCols)
                candidate = load(cfg.file, fileInfo(i).name);
                if isnumeric(candidate.(fileInfo(i).name))
                    variableName = fileInfo(i).name;
                    break;
                end
            end
        end
        if isempty(variableName)
            error('无法在 %s 中找到至少包含 %d 列的数值型 IMU 矩阵。', cfg.file, neededCols);
        end
    end

    loaded = load(cfg.file, variableName);
    data = loaded.(variableName);
    if ~isnumeric(data) || ndims(data) ~= 2
        error('MAT 变量 %s 必须是二维数值矩阵。', variableName);
    end
    if size(data,2) < max([cfg.accelCols(:); cfg.gyroCols(:)]) && size(data,1) >= max([cfg.accelCols(:); cfg.gyroCols(:)])
        data = data.';
    end
    if size(data,2) < max([cfg.accelCols(:); cfg.gyroCols(:)])
        error('MAT 变量 %s 不包含当前配置指定的 IMU 列。', variableName);
    end

    reader.data = data;
    reader.index = 1;
    reader.variableName = variableName;
    reader.cfg = cfg;
    fprintf('使用真实 IMU 的 MAT 变量 %s，共 %d 个采样点。\n', variableName, size(data,1));
end

function [wm, vm, isEof, reader] = readRealImuMatBlock(reader, sampleCount, rawTs, g0)
    startIndex = reader.index;
    endIndex = startIndex + sampleCount - 1;
    if endIndex > size(reader.data,1)
        wm = zeros(0,3);
        vm = zeros(0,3);
        isEof = true;
        return;
    end

    raw = double(reader.data(startIndex:endIndex,:));
    cfg = reader.cfg;
    accel = raw(:,cfg.accelCols) - cfg.accelBiasMps2(:).';
    if strcmpi(cfg.accelUnit, 'g')
        accel = accel*g0;
    elseif ~strcmpi(cfg.accelUnit, 'mps2')
        error('imuCfg.accelUnit 必须设置为 ''g'' 或 ''mps2''。');
    end
    gyroRadPerSec = raw(:,cfg.gyroCols) - cfg.gyroBiasRadPerSec(:).';

    accelPlatform = (cfg.CdataToPlatform*accel.').';
    gyroPlatform = (cfg.CdataToPlatform*gyroRadPerSec.').';
    vm = sum(accelPlatform*rawTs, 1);
    wm = sum(gyroPlatform*rawTs, 1);
    reader.index = endIndex + 1;
    isEof = false;
end

function reader = initRealImuReader(cfg)
    fid = fopen(cfg.file, 'r');
    if fid < 0
        error('无法打开真实 IMU 数据文件：%s', cfg.file);
    end
    for i = 1:cfg.headerLines
        if feof(fid)
            error('表头行数超过了真实 IMU 数据文件的实际长度。');
        end
        fgetl(fid);
    end
    reader.fid = fid;
    reader.cfg = cfg;
end

function [wm, vm, isEof] = readRealImuBlock(reader, nn, ts, g0, arcdeg)
    cfg = reader.cfg;
    wm = zeros(nn,3);
    vm = zeros(nn,3);
    isEof = false;
    row = 0;
    needCols = max([cfg.accelCols(:); cfg.gyroCols(:)]);

    while row < nn
        if feof(reader.fid)
            isEof = true;
            wm = wm(1:row,:);
            vm = vm(1:row,:);
            return;
        end

        line = fgetl(reader.fid);
        values = parseNumericTextLine(line);
        if numel(values) < needCols
            continue;
        end

        row = row + 1;
        accelPulse = values(cfg.accelCols);
        gyroPulse = values(cfg.gyroCols);

        accG = accelPulse(:)' ./ (cfg.accelScalePulsePerSecPerG(:)'*ts) + cfg.accelBiasG(:)';
        gyroDeg = gyroPulse(:)' ./ cfg.gyroScalePulsePerDeg(:)' + cfg.gyroBiasDeg(:)';

        vm(row,:) = accG*g0*ts;
        wm(row,:) = gyroDeg*arcdeg;
    end
    wm = sum(wm,1);
    vm = sum(vm,1);
end

function values = parseNumericTextLine(line)
    if ~ischar(line) && ~isstring(line)
        values = [];
        return;
    end
    text = char(line);
    text = strrep(text, ',', ' ');
    text = strrep(text, ';', ' ');
    text = strrep(text, sprintf('\t'), ' ');
    values = sscanf(text, '%f').';
end
