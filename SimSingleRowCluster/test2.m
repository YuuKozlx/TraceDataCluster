clc; clear;

% 定义体积的尺寸
n_cols = 256;   % x 方向的体素数量
n_rows = 256;   % y 方向的体素数量
n_slices = 1024; % z 方向的体素数量

% 定义球体的半径
radius = 3;

% 定义球体的中心坐标（可以不在 ISO 中心）
center_x = 50;  % 可以偏离 ISO 中心
center_y = 50;
center_z = 64;

pos = 30
% 创建体积几何体
vol_geom = astra_create_vol_geom(n_cols, n_rows, n_slices,pos-16,pos+16,pos-16,pos+16,0,128);

% 创建一个 3D 矩阵来存储球体数据
volume_data = zeros(n_cols, n_rows, n_slices, 'single');

% 填充球体数据：在体积矩阵中生成球体
for x = 1:n_cols
    for y = 1:n_rows
        for z = 1:n_slices
            distance = sqrt((x - center_x)^2 + (y - center_y)^2 + (z - center_z)^2);
            if distance <= radius
                volume_data(x, y, z) = 1;
            else
                volume_data(x, y, z) = 0.1;
            end
        end
    end
end

% 创建 ASTRA 的 3D 体数据对象，并使用 volume_data 进行初始化
vol_id = astra_mex_data3d('create', '-vol', vol_geom, volume_data);

% 定义旋转角度
angles = linspace(-pi/4, 7*pi/4, 361);  % 180 个角度
det_row_count = 2176;   % 探测器行数
det_col_count = 2176;   % 探测器列数
source_iso_distance = 400;  % 源到 ISO 的固定距离
iso_det_distance = 300;  % ISO到探测器的距离

% 初始化一个 3D 矩阵来存储所有角度的投影数据
proj_data_all = zeros(det_row_count, det_col_count, length(angles), 'single');

% ISO 中心在体积坐标系中的位置
iso_center_x = n_cols / 2;
iso_center_y = n_rows / 2;
iso_center_z = n_slices / 2;

for i = 1:length(angles)
    % 计算每个角度的光源位置
    theta = angles(i);

    % 创建 3D 锥形束投影几何体
    proj_geom = astra_create_proj_geom('cone', 0.089, 0.089, det_row_count, det_col_count, theta,source_iso_distance , iso_det_distance);

    % 创建投影数据对象
    proj_id = astra_mex_data3d('create', '-proj3d', proj_geom);

    % 设置投影算法
    cfg = astra_struct('FP3D_CUDA');
    cfg.VolumeDataId = vol_id;
    cfg.ProjectionDataId = proj_id;

    % 运行投影
    fp_id = astra_mex_algorithm('create', cfg);
    astra_mex_algorithm('run', fp_id);

    % 提取当前角度的投影数据
    proj_data_all(:, :, i) = astra_mex_data3d('get', proj_id);

    % 清理数据对象和算法
    astra_mex_algorithm('delete', fp_id);
    astra_mex_data3d('delete', proj_id);
end

figure(1);
xslice = center_x;  % 在 x 方向切片
yslice = center_y;  % 在 y 方向切片
zslice = center_z;  % 在 z 方向切片
slice(volume_data, xslice, yslice, zslice);  % 绘制切片
shading interp
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Ball Visualization using Slice');


% 另一种可视化：使用 isosurface 函数绘制球体表面
figure(2);
p = patch(isosurface(volume_data, 0.5));  % 设置等值面
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
isonormals(volume_data, p);  % 设置法线以改善渲染效果
view(3);  % 设置为三维视角
axis vis3d; % 保持纵横比例
camlight;  % 添加光照效果
lighting phong;  % 使用 phong 光照模型
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Ball Visualization using Isosurface');


% 显示某个角度的投影图像
m = max(proj_data_all,[],[1,2]);
m = squeeze(m);

figure(3);
for i =1:3:361
subplot(2,2,1)
imagesc(proj_data_all(:, :, i)');  % 显示第90个角度的投影
colormap('gray');
colorbar;
xlabel('Detector X');
ylabel('Detector Y');
title("X-ray Projection of a Rotating 3D Sphere, i = " + num2str(i) );
axis equal;
pause(0.2);
subplot(2,2,2)
hold on
temp = proj_data_all(:, :, i);

[row,col] = find(temp == m(i));
scatter(row,col);
end

% 清理体数据对象
astra_mex_data3d('delete', vol_id);
