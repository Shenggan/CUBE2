# 快速入门
## 设置参数
### 在[parameters.f90](./parameters.f90)中设置基础参数
- 环境变量
    - XFLAG_NO_OMP=-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large
    - XFLAG=-O3 -fpp -qopenmp -coarray=distributed -mcmodel=large
    - FFTFLAG=-I/opt/intel/oneapi/mkl/latest/include/fftw/ -qmkl
- 模拟参数
    - opath: 
        模拟数据的目录，所有模拟数据快照将存在这里
    - nn: 
        每一个纬度的节点数
    - box: 
        模拟的盒子的大小/Mpc
    - ncore: 
        每一个节点的核心数 > ntram\*nnest
    - nteam: 
        第一层OpenMP并行数
    - nnest: 
        第二层OpenMP并行数
    - ng: 
        模拟的标准格点大小
    - nic: 
        初始条件的格点大小
    - z_checkpoint: 
        快照点的红移，读取自根目录的 [z_checkpoint.txt](./z_checkpoint.txt)
- 宇宙学参数
    - h0: 哈勃常数
    - Mass_nu: 中微子质量和
    - omega_cdm: 冷暗物质质量分数
    - omega_bar: 重子质量分数
    - omega_mhd: 热暗物质质量分数
### 在[initialize.f90](./initialize.f90)中设置
- read_Gks: 
    是否读取已有的格林函数，首次运行需为.false.
- sim%cur_checkpoint: 
    开始模拟的快照点
### 在[kick.f90](./kick.f90)中设置
- PM3: 是否计算PM3的力
- PP: 是否计算PP力
## 配置环境
- 在319服务器中需先运行[module_load_t7920.sh](./module_load_t7920.sh)
## 生产初始条件
- 进入utilities目录
```batch
cd utilities
```
- 编译初始条件程序ic.x
```batch
make ic.x
```
- 回到主目录运行ic.x或者提交ic.sh
```batch
cd ..
./utilities/ic.x
    or 
sbatch ic.sh
```
## 开始模拟
- 首次运行
    - 编译主程序main.x，并运行或提交main.sh
    ```batch
    make
    ./main.x
    or
    sbatch main.sh
    ```
- 断点运行
    - 根据输出或者数据目录确定已完成的快照z = z_checkpoint(i)
    - 设置[initialize.f90](./initialize.f90)  line:34
        ```fortran
        logical,parameter :: read_Gks=.true.

        ···

        sim%cur_checkpoint= i 
        ```
    - 编译运行
## 故障排查
## 数据处理
### 功率谱
- 设置需要的产物[cicpower.f90](./utilities/cicpower.f90)
    - 计算快照点z_checkpoint(i)到z_checkpoint(j)的功率谱
        ```fortran
        do cur_checkpoint= i,j
        ```
    - 使用gadget的模拟结构
        ``` fortran
        #define gadget
        ···
        ! read Gadget output
        open(11,file='path_to_gadget_output_file',access='stream')

        ```
    - 获取投影
        ```fortran
        #define merge_projection
        ```
    - 获取中微子密度场
        ```fortran
        #define density_nu
        ```
    - 获取物质密度场
        ```fortran
        #define density_matter
        ```
- 进行功率谱修正[powerspectrum.f90](./utilities/powerspectrum.f90)   解除注释
    ```fortran
    subroutine auto_power(xi,cube1,n_particle,n_interp)

    ···
    
    xi(5,:)=xi(4,:)
    call pk_correction(xi,n_interp,3)
    call pk_correction(xi,n_interp,3)
    call pk_correction(xi,n_interp,3)
    call pk_correction(xi,n_interp,3)
    ! divide and normalize
    
    ···

    endsubroutine auto_power
    ```
- 编译并运行cicpower.x（参考[ic.x](#生产初始条件)）
### 查找halo
- 编译并运行fof.x（参考[ic.x](#生产初始条件)）
