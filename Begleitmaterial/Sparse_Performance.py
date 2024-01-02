import numpy as np
import scipy as sp

from timeit import Timer
import matplotlib.pyplot as plt


""" Laufzeitvergleich für die Multiplikation dünnbesetzter Matrizen """

def Sparse_Performance(
    N: int # Maximale Problemgröße
):

	# Initialisierung
	dims  = range(1,N+1)
	times = np.zeros((N,5))

	# Iterationsschleife
	for i in dims:

		# Definition der Funktionen
		with_csr       = lambda: sp.sparse.identity(i,format='csr') @ np.ones(i)
		with_csc       = lambda: sp.sparse.identity(i,format='csc') @ np.ones(i)
		with_coo       = lambda: sp.sparse.identity(i,format='coo') @ np.ones(i)
		without_sparse = lambda: np.identity(i) @ np.ones(i)
		no_zeros       = lambda: np.ones(i) * np.ones(i)
		
		# Laufzeitmessung
		times[i-1,0] = Timer(with_csr).timeit(number=1)
		times[i-1,1] = Timer(with_csc).timeit(number=1)
		times[i-1,2] = Timer(with_coo).timeit(number=1)
		times[i-1,3] = Timer(without_sparse).timeit(number=1)
		times[i-1,4] = Timer(no_zeros).timeit(number=1)

	# Plot
	plt.figure(figsize=(6,3))
	plt.loglog(dims,times[:,0])
	plt.loglog(dims,times[:,1])
	plt.loglog(dims,times[:,2])
	plt.loglog(dims,times[:,3])
	plt.loglog(dims,times[:,4])
	plt.xlabel('$N$')
	plt.ylabel('Laufzeit in Sekunden')
	plt.legend([
		r'$\mathbf{I}_{N \times N}\cdot\mathbf{1}_{N \times 1}$ (mit CSR)',
		r'$\mathbf{I}_{N \times N}\cdot\mathbf{1}_{N \times 1}$ (mit CSC)',
		r'$\mathbf{I}_{N \times N}\cdot\mathbf{1}_{N \times 1}$ (mit COO)',
		r'$\mathbf{I}_{N \times N}\cdot\mathbf{1}_{N \times 1}$ (ohne Kompression)',
		r'$\mathbf{1}_{N \times 1}\odot\mathbf{1}_{N \times 1}$',
	])
	plt.grid(visible=True,which='both',axis='both',color='gray',alpha=.25)
	plt.tight_layout()
	plt.show()
