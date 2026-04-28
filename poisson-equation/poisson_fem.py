import matplotlib.pyplot as plt
import numpy as np

# Triangular domain
p1 = [0, 0]
p2 = [1, 0]
p3 = [1, 1]
plt.plot(*zip(p1, p2, p3, p1), color="black")

# Triangular mesh
plt.plot(*zip([0, 0], [1/4, 0], [1/4, 1/4], [0, 0]), color="black")
plt.plot(*zip([1/4, 0], [2/4, 0], [2/4, 1/4], [1/4, 0]), color="black")
plt.plot(*zip([2/4, 0], [3/4, 0], [3/4, 1/4], [2/4, 0]), color="black")
plt.plot(*zip([3/4, 0], [1, 0], [1, 1/4], [3/4, 0]), color="black")
plt.plot(*zip([1/4, 1/4], [2/4, 1/4], [2/4, 2/4], [1/4, 1/4]), color="black")
plt.plot(*zip([2/4, 1/4], [3/4, 1/4], [3/4, 2/4], [2/4, 1/4]), color="black")
plt.plot(*zip([3/4, 1/4], [1, 1/4], [1, 2/4], [3/4, 1/4]), color="black")
plt.plot(*zip([2/4, 2/4], [3/4, 2/4], [3/4, 3/4], [2/4, 2/4]), color="black")
plt.plot(*zip([3/4, 2/4], [1, 2/4], [1, 3/4], [3/4, 2/4]), color="black")
plt.plot(*zip([3/4, 3/4], [1, 3/4], [1, 1], [3/4, 3/4]), color="black")
plt.show()

# Global stiffness matrix assembly
nodos = 15
n = 3
matrizcero = np.zeros((nodos, nodos))
matrizrigidez = np.array([[0.5, -0.5, 0], [-0.5, 1, -0.5], [0, -0.5, 0.5]])

for i in range(n):
    for j in range(n):
        matrizcero[j][i] = matrizrigidez[j][i]


def e(nodoi, nodoj, nodok, nodos):
    global e1
    e1 = np.zeros((nodos, nodos))
    matrizrigidez = np.array([[0.5, -0.5, 0], [-0.5, 1, -0.5], [0, -0.5, 0.5]])
    e1[nodoi][nodoi] = matrizrigidez[0][0]
    e1[nodoi][nodoj] = matrizrigidez[0][1]
    e1[nodoi][nodok] = matrizrigidez[0][2]
    e1[nodoj][nodoi] = matrizrigidez[1][0]
    e1[nodoj][nodoj] = matrizrigidez[1][1]
    e1[nodoj][nodok] = matrizrigidez[1][2]
    e1[nodok][nodoi] = matrizrigidez[2][0]
    e1[nodok][nodoj] = matrizrigidez[2][1]
    e1[nodok][nodok] = matrizrigidez[2][2]
    return e1


E1  = e(4, 2, 1, 15)
E2  = e(1, 3, 4, 15)
E3  = e(7, 4, 3, 15)
E4  = e(3, 6, 7, 15)
E5  = e(11, 7, 6, 15)
E6  = e(6, 10, 11, 15)
E7  = e(2, 4, 5, 15)
E8  = e(8, 5, 4, 15)
E9  = e(4, 7, 8, 15)
E10 = e(12, 8, 7, 15)
E11 = e(7, 11, 12, 15)
E12 = e(5, 8, 9, 15)
E13 = e(13, 9, 8, 15)
E14 = e(8, 12, 13, 15)
E15 = e(9, 13, 14, 15)

ensamble = E1+E2+E3+E4+E5+E6+E7+E8+E9+E10+E11+E12+E13+E14+E15+matrizcero

K = np.zeros((10, 10))
d = 10
for i in range(d):
    for j in range(d):
        K[j][i] = ensamble[i][j]

# Linear system KU = f
A = K
ELS = 16
B = np.array([[1/(6*ELS)], [3/(6*ELS)], [3/(6*ELS)], [3/(6*ELS)],
              [6/(6*ELS)], [3/(6*ELS)], [3/(6*ELS)], [6/(6*ELS)],
              [6/(6*ELS)], [3/(6*ELS)]])

A = np.array(A, dtype=float)
AB = np.concatenate((A, B), axis=1)
AB0 = np.copy(AB)

tamano = np.shape(AB)
n = tamano[0]
m = tamano[1]

# Partial pivoting
for i in range(0, n-1, 1):
    columna = abs(AB[i:, i])
    dondemax = np.argmax(columna)
    if dondemax != 0:
        temporal = np.copy(AB[i, :])
        AB[i, :] = AB[dondemax+i, :]
        AB[dondemax+i, :] = temporal

# Forward elimination
for i in range(0, n-1, 1):
    pivote = AB[i, i]
    adelante = i + 1
    for k in range(adelante, n, 1):
        factor = AB[k, i] / pivote
        AB[k, :] = AB[k, :] - AB[i, :] * factor

# Backward substitution
ultfila = n - 1
ultcolumna = m - 1
for i in range(ultfila, 0-1, -1):
    pivote = AB[i, i]
    atras = i - 1
    for k in range(atras, 0-1, -1):
        factor = AB[k, i] / pivote
        AB[k, :] = AB[k, :] - AB[i, :] * factor
    AB[i, :] = AB[i, :] / AB[i, i]

X = np.copy(AB[:, ultcolumna])
X = np.transpose([X])
print('Augmented matrix:')
print(AB0)
print('Solution X:')
print(X)

# Solution visualization
z = np.array([[0, 0, 0, 0, 0],
              [0.13924632, 0.13265931, 0.11167279, 0.0714614, 0],
              [0.22916667, 0.21721814, 0.18007047, 0.11167279, 0],
              [0.28048407, 0.2644761,  0.21721814, 0.13265931, 0],
              [0.3013174,  0.28048407, 0.22916667, 0.13924632, 0]])

fig, ax = plt.subplots()
ax.imshow(z)
plt.show()
