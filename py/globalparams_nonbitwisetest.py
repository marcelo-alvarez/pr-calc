import pandas as pd
import matplotlib.pyplot as plt

#act = pd.read_csv('/Users/margaretikape/Desktop/QUALS/Thesis/code/cmbonly_spectra_dr4.01/act_dr4.01_D_ell_TT_cmbonly.txt', sep='\s+', header=None)
#act.columns = ['ells', 'Dells', 'sigma_Dells']
#plt.errorbar(act['ells'], ((act['Dells']-act['Dells']) /act['Dells'] ), yerr=act['sigma_Dells'], fmt='.', color='black', ecolor='black', elinewidth=3, capsize=0);

GlobalParams = pd.read_csv('/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/newtest.ksz_cl.txt', sep='\s+', header=None)
GlobalParams.columns = ['ell', 'cll', 'dell']

NoGlobalParams = pd.read_csv('/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc-check/c++/example/test.ksz_cl.txt', sep='\s+', header=None)
NoGlobalParams.columns = ['ell', 'cll', 'dell']

plt.plot(GlobalParams['ell'], ((GlobalParams['dell'] - NoGlobalParams['dell']) / NoGlobalParams['dell']), 'r', label='GlobalParams')
plt.plot(NoGlobalParams['ell'], ((NoGlobalParams['dell'] - NoGlobalParams['dell'])/NoGlobalParams['dell']), 'k', label='NoGlobalParams')

plt.ylim(-0.001, 0.001)
plt.xlabel(r'$\ell$', fontsize=20)
plt.ylabel(r'$\Delta \ell/\ell$', fontsize=20)

plt.legend()
plt.show()

plt.plot(GlobalParams['ell'], GlobalParams['dell'], 'r', label='GlobalParams')
plt.plot(NoGlobalParams['ell'], NoGlobalParams['dell'], '--k', label='NoGlobalParams')

plt.xlabel(r'$\ell$', fontsize=20)
plt.ylabel(r'$\Delta \ell/\ell$', fontsize=20)

plt.legend()
plt.show()

GlobalParams = pd.read_csv('/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/output/newtest.history', sep='\s+', header=None)
GlobalParams.columns = ['z', 'chi', 'fcoll']

NoGlobalParams = pd.read_csv('/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc-check/c++/example/output/test.history', sep='\s+', header=None)
NoGlobalParams.columns = ['z', 'chi', 'fcoll']

plt.plot(NoGlobalParams['z'], ((NoGlobalParams['chi'] - GlobalParams['chi']) / GlobalParams['chi']), 'r', label='NoGlobalParams')
plt.plot(GlobalParams['z'], ((GlobalParams['chi'] - GlobalParams['chi'])/GlobalParams['chi']), 'k', label='GlobalParams')

#plt.ylim(-0.1, 0.1)
plt.xlabel(r'z', fontsize=20)
plt.ylabel(r'$\Delta \chi/\chi$', fontsize=20)
plt.legend()
plt.show()

plt.plot(GlobalParams['z'], GlobalParams['chi'] , 'r', label='GlobalParams')
plt.plot(NoGlobalParams['z'], NoGlobalParams['chi'], 'k', label='NoGlobalParams')

plt.xlim(5, 20)
plt.xlabel('z', fontsize=20)
plt.ylabel(r'$\chi$', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend()
plt.show()

plt.plot(GlobalParams['z'], GlobalParams['fcoll'] , 'r', label='GlobalParams')
plt.plot(NoGlobalParams['z'], NoGlobalParams['fcoll'], 'k', label='NoGlobalParams')

plt.xlim(-0.1, 15)
plt.xlabel('z', fontsize=20)
plt.ylabel(r'$f_{\rm coll}$', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend()
plt.show()

plt.plot(GlobalParams['z'], ((GlobalParams['fcoll'] - NoGlobalParams['fcoll']) / NoGlobalParams['fcoll']), 'r', label='GlobalParams')
plt.plot(NoGlobalParams['z'], ((NoGlobalParams['fcoll'] - NoGlobalParams['fcoll'])/NoGlobalParams['fcoll']), 'k', label='NoGlobalParams')

plt.ylim(-0.00000001, 0.00000001)
plt.xlabel('z', fontsize=20)
plt.ylabel(r'$\Delta f_{coll}/f_{coll}$', fontsize=20)
plt.legend()
plt.show()

import numpy as np
import matplotlib.pyplot as plt

GlobalParams = open("/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/maps/newtest.kszmap","rb")
NoGlobalParams = open("/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc-check/c++/example/maps/test.kszmap","rb")

# Read size and Field of view
n = np.fromfile(GlobalParams, dtype=np.int32, count=2)
n2 = np.fromfile(NoGlobalParams, dtype=np.int32, count=2)

# Read image, need to convert to float64 (double) to avoid precision issues (affects mean and var)
GlobalParams = np.fromfile(GlobalParams, dtype=np.float32, count=n[0]*n[1]).astype(np.float64).reshape(n)
GlobalParams = GlobalParams[:n[0],:n[1]]
NoGlobalParams = np.fromfile(NoGlobalParams, dtype=np.float32, count=n2[0]*n2[1]).astype(np.float64).reshape(n2)
NoGlobalParams = NoGlobalParams[:n2[0],:n2[1]]

plt.plot(GlobalParams[:n[0]], (GlobalParams[:n[1]]-GlobalParams[:n[1]])/GlobalParams[:n[1]], 'or')#, label='Colossus')
plt.plot(NoGlobalParams[:n2[0]], (NoGlobalParams[:n2[1]]-GlobalParams[:n[1]])/GlobalParams[:n[1]], 'k')#, label='WMAP')
#plt.ylim(-0.00000001, 0.00000001)
plt.legend()
plt.show()
