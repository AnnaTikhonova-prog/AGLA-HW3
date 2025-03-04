
def gaussian_elimination(matrix):
    n = len(matrix)
    m = len(matrix[0])
    matrix = [row[:] for row in matrix]
    lead = 0

    for r in range(n):
        if lead >= m:
            return matrix, lead
        i = r
        while matrix[i][lead] == 0:
            i += 1
            if i == n:
                i = r
                lead += 1
                if lead == m:
                    return matrix, lead

        matrix[i], matrix[r] = matrix[r], matrix[i]
        lv = matrix[r][lead]
        for j in range(n):
            if j != r:
                factor = matrix[j][lead] / lv
                for k in range(lead, m):
                    matrix[j][k] -= factor * matrix[r][k]
        lead += 1
    return matrix, lead


def rank(matrix):
    ref, _ = gaussian_elimination(matrix)
    return sum(1 for row in ref if any(row))


def nullspace_basis(matrix):
    n = len(matrix)
    m = len(matrix[0])
    ref, lead = gaussian_elimination(matrix)
    rank_A = rank(matrix)

    if rank_A == m:
        return []

    pivot_cols = []
    for r in range(rank_A):
        for i in range(m):
            if ref[r][i] != 0:
                pivot_cols.append(i)
                break

    free_vars = [i for i in range(m) if i not in pivot_cols]
    basis = []

    for free in free_vars:
        vector = [0] * m
        vector[free] = 1
        for r in range(rank_A):
            pivot = next(i for i in range(m) if ref[r][i] != 0)
            if free > pivot:
                vector[pivot] = -ref[r][free] / ref[r][pivot]
        basis.append(vector)
    return basis

def transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]


def find_x(matrix):
    basis = nullspace_basis(matrix)
    if not basis:
        return [0] * len(matrix[0])
    return basis[0]

def find_y(matrix):
    A_T = transpose(matrix)
    basis = nullspace_basis(A_T)
    if not basis:
        return [0] * len(A_T[0])
    return basis[0]


def find_z(matrix):
    for row in matrix:
        if any(row):
            return row
    return [0] * len(matrix[0])



A = [[1, 2, 1],
     [2, 4, 3],
     [3, 6, 4]]

x = find_x(A)
y = find_y(A)
z = find_z(A)

print("x orthogonal to row space:", x)
print("y orthogonal to coolumn space:", y)
print("z orthogonal to null space:", z)