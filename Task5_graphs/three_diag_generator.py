from scipy.sparse import diags
import numpy
from scipy import linalg

n = 100
rand_diag = numpy.random.randint(100, 200, n)
rand_subdiag = numpy.random.randint(1, 30, n - 1)
rand_subdiag1 = numpy.random.randint(1, 30, n - 2)
rand_subdiag2 = numpy.random.randint(1, 30, n - 3)

k = [numpy.ones(n - 3) * rand_subdiag2, numpy.ones(n - 2) * rand_subdiag1,
     numpy.ones(n - 1) * rand_subdiag, numpy.ones(n) * rand_diag, numpy.ones(n - 1) * rand_subdiag,
     numpy.ones(n - 2) * rand_subdiag1, numpy.ones(n - 3) * rand_subdiag2]
offset = [-3, -2, -1, 0, 1, 2, 3]
A = diags(k, offset).toarray()
numpy.savetxt("matrix.txt", A)
eigs = linalg.eigh(A, eigvals_only = True)
to_write = numpy.array((numpy.min(eigs), numpy.max(eigs)), dtype="str")

with open("matrix.txt", "a") as file:
    file.write(" ".join(to_write))
