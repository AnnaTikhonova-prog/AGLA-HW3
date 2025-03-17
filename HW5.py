import math

def dot_product(v1, v2):
    return sum(x * y for x, y in zip(v1, v2))

def vector_subtract(v1, v2):
    return [x - y for x, y in zip(v1, v2)]

def scalar_multiply(scalar, vector):
    return [scalar * x for x in vector]

def normalize(vector):
    norm = math.sqrt(dot_product(vector, vector))
    return [x / norm for x in vector]

def gram_schmidt(vectors):
    basis = []
    for v in vectors:
        w = v[:]
        for b in basis:
            projection = scalar_multiply(dot_product(w, b) / dot_product(b, b), b)
            w = vector_subtract(w, projection)
        w_norm = math.sqrt(dot_product(w, w))
        if w_norm > 1e-10:
            w_normalized = normalize(w)
            basis.append(w_normalized)
    return basis

vectors = [
    [1, 2, 3, 4],
    [2, 3, 4, 5],
    [7, 8, 5, 6],
    [9, 2, 6, 7]
]

orthonormal_basis = gram_schmidt(vectors)
for i, vector in enumerate(orthonormal_basis):
    print(f"Vector {i+1}: {vector}")