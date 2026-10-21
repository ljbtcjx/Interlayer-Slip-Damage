import os

input_file = "relaxed_initial.xyz"
output_step0 = "step_0_standard.xyz" # 修正格式后的初始态 (备用)
output_step1 = "step_1.xyz"          # 拉伸后的第一步 (用来跑任务)

# 1. 设置拉伸参数
# 根据你的数据，Z坐标在21-26的顶部碳层是 145~184 号原子
pull_indices = range(145, 184 + 1)  
delta_z = 0.2  # 第一步拉伸步长 0.2 Å

with open(input_file, 'r') as f:
    lines = f.readlines()

num_atoms = lines[0].strip()
comment = "Standardized XYZ for CP2K - Cf/SiC Interface"

def write_xyz(filename, lines, is_pulled=False):
    with open(filename, 'w') as f:
        f.write(f"{num_atoms}\n")
        f.write(f"{comment}\n")
        for i, line in enumerate(lines[2:]):
            parts = line.split()
            if len(parts) >= 4:
                # 提取 OVITO 格式的坐标和元素
                x, y, z, element = float(parts[0]), float(parts[1]), float(parts[2]), parts[3]
                
                atom_id = i + 1 
                
                # 如果是拉伸文件，并且当前原子属于顶层，则 Z 坐标增加
                if is_pulled and atom_id in pull_indices:
                    z += delta_z
                    
                # 按照 CP2K 要求的标准格式写入：元素 X Y Z
                f.write(f"{element:3s} {x:15.10f} {y:15.10f} {z:15.10f}\n")

# 生成标准格式初始态
write_xyz(output_step0, lines, is_pulled=False)
# 生成拉伸第一步
write_xyz(output_step1, lines, is_pulled=True)

print("🎉 格式修正与拉伸操作完成！")
print(f"👉 标准化初始态已生成: {output_step0}")
print(f"👉 拉伸 0.2 Å 的计算文件已生成: {output_step1}")