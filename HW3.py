
def round_matrix(matrix):
    return [[round(val, 2) for val in row] for row in matrix]



def input_matrix():
    rows = int(input("Enter number of rows: "))
    cols = int(input("Enter number of columns "))
    matrix = []

    for i in range(rows):
        row = list(map(float, input(f"Enter a row {i + 1} (elements separated by space): ").split()))
        if len(row) != cols:
            raise ValueError("The number of elements in the row doesnt match the number of columns.")
        matrix.append(row)
    return matrix


def reduced_row_echelon(matrix, epsilon=1e-12):
    mat = [row.copy() for row in matrix]
    row_count, col_count = len(mat), len(mat[0])
    pivots = []
    lead = 0
    for r in range(row_count):
        if lead >= col_count:
            break
        pivot = max(range(r, row_count), key=lambda i: abs(mat[i][lead]))
        if abs(mat[pivot][lead]) < epsilon:
            lead += 1
            continue
        mat[r], mat[pivot] = mat[pivot], mat[r]
        pivot_val = mat[r][lead]
        mat[r] = [val / pivot_val for val in mat[r]]
        for i in range(row_count):
            if i != r:
                factor = mat[i][lead]
                mat[i] = [a - factor * b for a, b in zip(mat[i], mat[r])]
        pivots.append(lead)
        lead += 1
    return mat, pivots


def extract_row_space(rref_mat):
    return [row for row in rref_mat if any(abs(val) > 1e-12 for val in row)]


def extract_col_space(original, pivots):
    return [[row[col] for row in original] for col in pivots]


def extract_null_space(rref_mat, pivots):
    cols = len(rref_mat[0])
    free_vars = [col for col in range(cols) if col not in pivots]
    null_basis = []

    for free_var in free_vars:
        vector = [0] * cols
        vector[free_var] = 1
        for r, pivot_col in enumerate(pivots):
            vector[pivot_col] = -sum(rref_mat[r][j] * vector[j] for j in range(cols) if j != pivot_col)
        null_basis.append(vector)
    return null_basis


def transpose(mat):
    return [list(row) for row in zip(*mat)]


def find_all_subspaces(matrix):
    rref_mat, pivots = reduced_row_echelon(matrix)
    col_space = extract_col_space(matrix, pivots)
    row_space = extract_row_space(rref_mat)
    null_space = extract_null_space(rref_mat, pivots)
    transposed = transpose(matrix)
    rref_t, pivots_t = reduced_row_echelon(transposed)
    row_null_space = extract_null_space(rref_t, pivots_t)
    return col_space, row_space, null_space, row_null_space


if __name__ == "__main__":
    matrix = input_matrix()
    c_space, r_space, cn_space, rn_space = find_all_subspaces(matrix)
    c_space = round_matrix(c_space)
    r_space = round_matrix(r_space)
    cn_space = round_matrix(cn_space)
    rn_space = round_matrix(rn_space)

    print("Column Space:", c_space)
    print("Row Space:", r_space)
    print("Column Null Space:", cn_space)
    print("Row Null Space:", rn_space)


# Example of work:
#
# Введите количество строк: 3
# Введите количество столбцов 3
# Введите строку 1 (элементы через пробел): 1 2 3
# Введите строку 2 (элементы через пробел): 4 5 6
# Введите строку 3 (элементы через пробел): 7 8 9
# Column Space: [[1.0, 4.0, 7.0], [2.0, 5.0, 8.0]]
# Row Space: [[1.0, 0.0, -1.0], [0.0, 1.0, 2.0]]
# Column Null Space: [[1.0, -2.0, 1]]
# Row Null Space: [[1.0, -2.0, 1]]

# Description of all functions:
# round_matrix(): Rounds all elements of a matrix to two decimal places.
# input_matrix(): Prompt the user to input a matrix row by row .
# reduced_row_echelon(): Convert the matrix to reduced row echelon form and records pivot columns.
# extract_row_space(): Extracts the row space.
# extract_col_space(): Extracts the column space using pivot columns from the original matrix.
# extract_null_space(): Finds the null space by solving for free variables.
# transpose(): Returns the transpose  matrix.
# find_all_subspaces(): Collects the column space, row space, null space, and left null space.
#