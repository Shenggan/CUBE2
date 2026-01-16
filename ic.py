import camb,sys,os,re
import numpy as np

def match_para(para):
    # 读取Fortran文件
    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()

    # 使用正则表达式匹配模式
    pattern = r'parameter\s*::\s*%s\s*=\s*([^\s]+)'%para
    match = re.search(pattern, content)

    if match:
        variable_value = match.group(1).strip()
        return float(variable_value)
    else:
        print(f"Pattern({para}) not found in the file.")
        sys.exit()

#get cube parm
if (1):
    H0=match_para('h0')*100
    omega_bar=match_para('omega_bar')
    omega_cdm=match_para('omega_cdm')
    omk=0.0

    ns=match_para('n_s')
    As=match_para('A_s')
    z_list = np.loadtxt('./z_checkpoint.txt')
    z_max = z_list[0]

    ombh2 = omega_bar*(H0/100)**2
    omch2 = omega_cdm*(H0/100)**2

    file_path = './parameters.f90'
    with open(file_path, 'r') as file:
        content = file.read()

    # 使用正则表达式匹配模式
    pattern = r'parameter\s*::\s*opath\s*=\s*([^\s]+)'
    match = re.search(pattern, content)
    opath = os.path.expanduser(match.group(1).strip()[1:-1])
    print(f"Variable value of opath : {opath}")
    print('\n'+('+'*40+'\n')*2)
    print('Cosmology  Paras:\n\n   omega_b:   %.6f\n   omega_c:   %.6f\n       A_s:   %.2e\n       n_s:   %.6f\n'%(omega_bar,omega_cdm,As,ns))



#mkdir
try:
    os.system('mkdir -p '+opath+'/IC')
except:
    None
    

npbin = 219  
k_ic_min = 1e-4
k_ic_max = 1e2

par = camb.CAMBparams()
par.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk)
par.InitPower.set_params(As=As,ns=ns)
par.set_matter_power(redshifts=z_list, kmax=k_ic_max, nonlinear=True)
par.NonLinear = camb.model.NonLinear_both
result = camb.get_results(par)
kh_ic, z_list,Pk_cb= result.get_matter_power_spectrum(minkh=k_ic_min, maxkh=k_ic_max, npoints = npbin,var1='delta_nonu',var2='delta_nonu')

for i  in range(len(z_list)):
    z = z_list[i]
    print('write %sIC/Pcb_%3.4f.txt'%(opath,z))
    np.savetxt(opath+f'/IC/Pcb_{z:.4f}.txt',np.array([kh_ic,Pk_cb[i]]).T)

Pk_cb_ic = Pk_cb[-1]
print('Pcb_ic saved in %sIC/Pcb_ic.txt'%opath)
np.savetxt(opath+'/IC/Pcb_ic.txt',np.array([kh_ic,Pk_cb_ic]).T)

