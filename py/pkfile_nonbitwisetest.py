import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#power spectrum percentage differences
wmappkfile = pd.read_csv('/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/newtest.ksz_cl.txt', sep='\s+', header=None)
wmappkfile.columns = ['ell', 'cll', 'dell']

colpkfile = pd.read_csv('/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/test.ksz_cl.txt', sep='\s+', header=None)
colpkfile.columns = ['ell', 'cll', 'dell']

plt.plot(colpkfile['ell'], ((colpkfile['dell'] - wmappkfile['dell']) / wmappkfile['dell']), 'r', label='Colossus p(k)')
plt.plot(wmappkfile['ell'], ((wmappkfile['dell'] - wmappkfile['dell'])/wmappkfile['dell']), 'k', label='wmap p(k)')

plt.ylim(-0.5, 0.5)
plt.xlabel(r'$\ell$', fontsize=20)
plt.ylabel(r'$\Delta \ell/\ell$', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend(loc='best', fontsize=20)
#plt.savefig('pkfile_nonbitwisetest_PSdiff.png',bbox_inches='tight')
plt.show()

#Power spectrum
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.semilogy(colpkfile['ell'], colpkfile['dell'], 'r', label='Colossus p(k)')
plt.semilogy(wmappkfile['ell'], wmappkfile['dell'], '--k', label='wmap p(k)')

plt.xlabel(r'$\ell$', fontsize=20)
plt.ylabel(r'$\Delta \ell/\ell$', fontsize=20)

plt.legend()
#plt.savefig('pkfile_nonbitwisetest_PS.png',bbox_inches='tight')
plt.show()

#chi percentage differences
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

wmappkfile = pd.read_csv('/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/output/newtest.history', sep='\s+', header=None)
wmappkfile.columns = ['z', 'chi', 'fcoll']

colpkfile = pd.read_csv('/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/output/test.history', sep='\s+', header=None)
colpkfile.columns = ['z', 'chi', 'fcoll']

plt.plot(colpkfile['z'], ((colpkfile['chi'] - wmappkfile['chi']) / wmappkfile['chi']), 'r', label='Colossus p(k)')
plt.plot(wmappkfile['z'], ((wmappkfile['chi'] - wmappkfile['chi'])/wmappkfile['chi']), 'k', label='wmap p(k)')

plt.xlabel(r'z', fontsize=20)
plt.ylabel(r'$\Delta \chi/\chi$', fontsize=20)
plt.legend()
#plt.savefig('pkfile_nonbitwisetest_CHIdiff.png',bbox_inches='tight')
plt.show()

#chi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.plot(colpkfile['z'], colpkfile['chi'] , 'r', label='Colossus p(k)')
plt.plot(wmappkfile['z'], wmappkfile['chi'], 'k', label='wmap p(k)')

plt.xlim(5, 20)
plt.xlabel('z', fontsize=20)
plt.ylabel(r'$\chi$', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend()
#plt.savefig('pkfile_nonbitwisetest_CHI.png',bbox_inches='tight')
plt.show()

#fcoll
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.semilogy(colpkfile['z'], colpkfile['fcoll'] , 'r', label='Colossus p(k)')
plt.semilogy(wmappkfile['z'], wmappkfile['fcoll'], 'k', label='wmap p(k)')

plt.xlim(-0.1, 15)
plt.xlabel('z', fontsize=20)
plt.ylabel(r'$f_{\rm coll}$', fontsize=20)
plt.tick_params(labelsize=20)
plt.legend()
#plt.savefig('pkfile_nonbitwisetest_fcoll.png',bbox_inches='tight')
plt.show()

#fcoll percenttage differences
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.plot(colpkfile['z'], ((colpkfile['fcoll'] - wmappkfile['fcoll']) / wmappkfile['fcoll']), 'r', label='Colossus p(k)')
plt.plot(wmappkfile['z'], ((wmappkfile['fcoll'] - wmappkfile['fcoll'])/wmappkfile['fcoll']), 'k', label='wmap p(k)')

plt.ylim(-0.00000001, 0.00000001)
plt.xlabel('z', fontsize=20)
plt.ylabel(r'$\Delta f_{coll}/f_{coll}$', fontsize=20)
plt.legend()
#plt.savefig('pkfile_nonbitwisetest_fcollDiff.png',bbox_inches='tight')
plt.show()

#T map diff
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mapfile = open("/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/maps/test.kszmap","rb")
mapfile2 = open("/global/cscratch1/sd/ikapem/ksz-reionization/pr-calc/c++/example/maps/newtest.kszmap","rb")

# Read size and Field of view
n = np.fromfile(mapfile, dtype=np.int32, count=2)
n2 = np.fromfile(mapfile2, dtype=np.int32, count=2)

# Read image, need to convert to float64 (double) to avoid precision issues (affects mean and var)
mapfile = np.fromfile(mapfile, dtype=np.float32, count=n[0]*n[1]).astype(np.float64).reshape(n)
mapfile = mapfile[:n[0],:n[1]]
mapfile2 = np.fromfile(mapfile2, dtype=np.float32, count=n2[0]*n2[1]).astype(np.float64).reshape(n2)
mapfile2 = mapfile[:n2[0],:n2[1]]

plt.plot(mapfile[:n[0]], (mapfile[:n[1]]-mapfile[:n[1]])/mapfile[:n[1]], 'or')#, label='Colossus')
plt.plot(mapfile2[:n2[0]], (mapfile2[:n2[1]]-mapfile[:n[1]])//mapfile[:n[1]], 'k')#, label='WMAP')
plt.ylim(-0.00000001, 0.00000001)
plt.legend()
#plt.savefig('pkfile_nonbitwisetest_PSmap.png',bbox_inches='tight')
plt.show()
