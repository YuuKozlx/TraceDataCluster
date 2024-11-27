% 定义体积的尺寸和球体半径
n_cols = 64;   % x 方向的体素数量
n_rows = 64;   % y 方向的体素数量
n_slices = 64; % z 方向的体素数量
radius = 20;   % 球体半径

% 创建体积几何体，使用 astra_create_vol_geom 创建 3D 几何体定义
vol_geom = astra_create_vol_geom(n_cols, n_rows, n_slices);

% 创建一个 3D 矩阵来存储球体数据
volume_data = zeros(n_cols, n_rows, n_slices, 'single');  % 确保数据类型为单精度

% 计算球体的中心坐标
center_x = n_cols / 2;
center_y = n_rows / 2;
center_z = n_slices / 2;

% 填充球体数据：在体积矩阵中生成球体
for x = 1:n_cols
    for y = 1:n_rows
        for z = 1:n_slices
            % 计算当前体素到中心的距离
            distance = sqrt((x - center_x)^2 + (y - center_y)^2 + (z - center_z)^2);
            
            % 如果在球体半径范围内，设置体素值为 1
            if distance <= radius
                volume_data(x, y, z) = 1;
            end
        end
    end
end

% 创建 ASTRA 的 3D 体数据对象，并使用 volume_data 进行初始化
vol_id = astra_mex_data3d('create', '-vol', vol_geom, volume_data);

% 可视化：使用 slice 函数在 x、y、z 方向上切片查看内部结构
figure;
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
figure;
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

% 清理数据对象
astra_mex_data3d('delete', vol_id);
