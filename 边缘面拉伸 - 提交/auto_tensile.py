import os
import time

# ==========================================
# 1. 全局核心参数设置 (根据你的需求修改)
# ==========================================
START_STEP = 12          # 从第几步开始跑 (你已经跑完了 Step 1)
END_STEP = 20           # 跑到第几步结束 (总共拉伸 20 步，即 4.0 Å)
DELTA_Z = 0.2           # 每步拉伸的步长 (Å)
BASE_CELL_Z = 34.4741712485839216  # 原始未拉伸前的晶胞 Z 轴长度

# 顶层被拉拽的碳原子编号 (145 到 154)
TOP_ATOMS = list(range(145, 155))

# 运行 CP2K 的命令 (请确保与你服务器的实际运行命令一致)
MPI_CMD = "mpirun -np 30 cp2k.popt"

# ==========================================
# 2. CP2K INP 文件模板
# ==========================================
def generate_inp_content(step):
    prev_step = step - 1
    current_cell_z = BASE_CELL_Z + step * DELTA_Z
    
    # 注意这里使用三引号定义多行字符串，大括号内是动态替换的变量
    inp_text = f"""&GLOBAL
  PROJECT Cf_SiC_Tensile_Step{step}
  RUN_TYPE GEO_OPT
  PRINT_LEVEL LOW
&END GLOBAL

&MOTION
  &GEO_OPT
    TYPE MINIMIZATION
    MAX_ITER 100
    OPTIMIZER BFGS
  &END GEO_OPT

  &CONSTRAINT
    &FIXED_ATOMS
      ! --- 底部固定的 SiC 基底 ---
      LIST 1 7 13 19 25 31 37 43 49 55 61 67 73 79 85
      LIST 91 97 103 109 115 121 127 133 139 185 186 191 192 197 198
      LIST 203 204 209 210 215 216 221 222 227 228 233 234 239 240 245
      LIST 246 251 252 257 258 263 264 269 270 275 276 281 282 287 288
      LIST 293 294 299 300 305 306 311 312 317 318 323 324
      ! --- 顶部被拉拽的碳层 ---
      LIST 145 146 147 148 149 150 151 152 153 154
    &END FIXED_ATOMS
  &END CONSTRAINT
&END MOTION

&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME POTENTIAL
    
    ! 🎯 动态继承上一步的波函数
    WFN_RESTART_FILE_NAME Cf_SiC_Tensile_Step{prev_step}-RESTART.wfn
    
    UKS .TRUE.
    MULTIPLICITY 1
    
    &MGRID
      CUTOFF 400
      REL_CUTOFF 60
    &END MGRID
    
    &SCF
      MAX_SCF 200
      EPS_SCF 1.0E-6
      SCF_GUESS RESTART
      ADDED_MOS 100
      &SMEAR
        METHOD FERMI_DIRAC
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.4
        NBROYDEN 8
      &END MIXING
      &OUTER_SCF
        MAX_SCF 10
        EPS_SCF 1.0E-6
      &END OUTER_SCF
    &END SCF
    
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
         POTENTIAL_TYPE PAIR_POTENTIAL
         &PAIR_POTENTIAL
            PARAMETER_FILE_NAME dftd3.dat
            TYPE DFTD3
            REFERENCE_FUNCTIONAL PBE
         &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT

  &SUBSYS
    &CELL
      A   12.3319422638933904    0.0000000000000000    0.0000000000000000
      B    0.0000000000000000   16.0196629178019876    0.0000000000000000
      ! 🎯 动态加长 C 轴
      C    0.0000000000000000    0.0000000000000000   {current_cell_z:.16f}
    &END CELL
    
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME ./step_{step}.xyz  
    &END TOPOLOGY
    
    &KIND C
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q4
    &END KIND
    
    &KIND Si
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
"""
    return inp_text

# ==========================================
# 3. 自动化主循环
# ==========================================
print(f"🚀 开启自动拉伸任务，从 Step {START_STEP} 运行至 Step {END_STEP}")

for step in range(START_STEP, END_STEP + 1):
    prev_step = step - 1
    prev_traj_file = f"Cf_SiC_Tensile_Step{prev_step}-pos-1.xyz"
    prev_wfn_file = f"Cf_SiC_Tensile_Step{prev_step}-RESTART.wfn"
    
    # 检查上一步是否成功跑完 (有没有生成轨迹文件和波函数)
    if not os.path.exists(prev_traj_file):
        print(f"❌ 致命错误：找不到上一步的轨迹文件 {prev_traj_file}。程序终止。")
        break
    if not os.path.exists(prev_wfn_file):
        print(f"⚠️ 警告：找不到上一步的波函数 {prev_wfn_file}，SCF 收敛可能会变慢。")

    print(f"==========================================")
    print(f"🛠️  正在准备 Step {step} 的计算...")
    
    # 3.1 提取上一步的最后一帧并拉伸顶层原子
    with open(prev_traj_file, 'r') as f:
        lines = f.readlines()
    
    natoms = int(lines[0].strip())
    # 截取整个轨迹文件的最后 (natoms + 2) 行，即最后一帧的基态
    last_frame = lines[-(natoms + 2):]
    
    new_xyz_filename = f"step_{step}.xyz"
    with open(new_xyz_filename, 'w') as f:
        f.write(last_frame[0])
        f.write(f"Cf_SiC interface stretched to step {step} (+{step * DELTA_Z} A)\n")
        
        for i in range(2, len(last_frame)):
            parts = last_frame[i].split()
            element = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            
            atom_id = i - 1
            if atom_id in TOP_ATOMS:
                z += DELTA_Z  # 向上拉长 0.2 A
                
            f.write(f"{element:3s} {x:15.10f} {y:15.10f} {z:15.10f}\n")
            
    print(f"✅ 成功提取并拉伸坐标，已生成 -> {new_xyz_filename}")
    
    # 3.2 生成新的 INP 输入文件
    new_inp_filename = f"step_{step}.inp"
    with open(new_inp_filename, 'w') as f:
        f.write(generate_inp_content(step))
    print(f"✅ 成功更新输入脚本，已生成 -> {new_inp_filename}")
    
    # 3.3 提交 CP2K 计算
    out_filename = f"step_{step}.out"
    full_cmd = f"{MPI_CMD} -i {new_inp_filename} > {out_filename}"
    print(f"⏳ 正在后台运行计算：{full_cmd} ...")
    
    # 调用系统命令执行 CP2K (程序会在这里阻塞，直到 CP2K 算完才继续往下走)
    os.system(full_cmd)
    
    # 3.4 提取并打印能量
    print(f"🎉 Step {step} 计算完成！")
    os.system(f"grep 'ENERGY| Total FORCE_EVAL' {out_filename} | tail -n 1")
    print("\n")

print("🏆 所有拉伸任务已全部自动化完成！快去画能量曲线吧！")