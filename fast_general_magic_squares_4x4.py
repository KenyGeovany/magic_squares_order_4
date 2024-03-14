import numpy as np
from datetime import datetime
from itertools import permutations

output_directory = 'output'

class magic_square:
    def __init__(self):
        # All elements must be different.
        self.elements = {0, 1, 2}
        # Default sum if the magic square is isotropic.
        self.s = 0
        # Customized sums by safes
        self.sR_0 = 0
        # [0,0] [0,1] [0,2] [0,3] \\
        #   *     *     *     *   \\
        #   *     *     *     *   \\
        #   *     *     *     *
        self.sR_1 = 0
        #   *     *     *     *   \\
        # [1,0] [1,1] [1,2] [1,3] \\
        #   *     *     *     *   \\
        #   *     *     *     *
        self.sR_2 = 0
        #   *     *     *     *   \\
        #   *     *     *     *   \\
        # [2,0] [2,1] [2,2] [2,3] \\
        #   *     *     *     *
        self.sC_0 = 0
        # [0,0]   *     *     *   \\
        # [1,0]   *     *     *   \\
        # [2,0]   *     *     *   \\
        # [3,0]   *     *     *
        self.sC_1 = 0
        #   *   [0,1]   *     *   \\
        #   *   [1,1]   *     *   \\
        #   *   [2,1]   *     *   \\
        #   *   [3,1]   *     *
        self.sC_2 = 0
        #   *     *   [0,2]   *   \\
        #   *     *   [1,2]   *   \\
        #   *     *   [2,2]   *   \\
        #   *     *   [3,2]   *
        self.sC_3 = 0
        #   *     *     *   [0,3] \\
        #   *     *     *   [1,3] \\
        #   *     *     *   [2,3] \\
        #   *     *     *   [3,3]
        self.sD_0 = 0
        # [0,0]   *     *     *   \\
        #   *   [1,1]   *     *   \\
        #   *     *   [2,2]   *   \\
        #   *     *     *   [3,3]
        self.sD_1 = 0
        #   *     *     *   [0,3] \\
        #   *     *   [1,2]   *   \\
        #   *   [2,1]   *     *   \\
        # [3,0]   *     *     *

        # Magic square type
        # If Magic square is not initialized then functions give error
        self.initialized = False
        # Isotropic Magic Square has each safe sum = s
        self.isotropic = False
        # Ordinary Magic Square has elements = {0,...,15}
        self.ordinary = False

        # Smallest element
        self.x0 = min(self.elements)

        # Biggest element
        self.y0 = max(self.elements)

    def set_isotropic(self):
        self.sR_0 = self.s
        self.sR_1 = self.s
        self.sR_2 = self.s
        self.sC_0 = self.s
        self.sC_1 = self.s
        self.sC_2 = self.s
        self.sC_3 = self.s
        self.sD_0 = self.s
        self.sD_1 = self.s
        self.x0 = min(self.elements)
        self.y0 = max(self.elements)
        self.isotropic = True

    def set_ordinary(self):
        return 0

    def compute_squares_of_type_I(self, show=False, save=False):
        # Magic square of type I
        #  x0  A[0] A[1] B[0] \\
        # A[2] A[3] A[4] B[1] \\
        # B[5] B[4] A[5] B[3] \\
        # B[6] B[7] B[8] B[2]
        if show:
            print("Computing squares of Type I \n", flush=True)
        # Remove x0 from the set of element possibilities
        domain = self.elements.difference({self.x0})
        # Compute all combinations of 6 elements taken from domain (15 elements)
        comb = permutations(domain, 6)
        # Output set. This set contains also no ordinary magic squares.
        squares = set()
        for A in comb:
            # Compute the complement of A u {x0} (9 elements)
            C = domain - set(A)
            # Compute the 9 dependent variables
            B = (self.sR_0 - (self.x0 + A[0] + A[1]),
                 self.sR_1 - (A[2] + A[3] + A[4]),
                 self.sD_0 - (self.x0 + A[3] + A[5]),
                 self.sC_3 - (self.sR_0 + self.sR_1 + self.sD_0) + (
                             2 * self.x0 + A[0] + A[1] + A[2] + 2 * A[3] + A[4] + A[5]),
                 int(0.5 * (self.sD_0 + self.sD_1 + self.sR_1 + self.sR_2 - self.sC_0 - self.sC_3)) - (
                             A[3] + A[4] + A[5]),
                 int(0.5 * (self.sD_0 - self.sD_1 + 2 * self.sR_0 + self.sR_1 + self.sR_2 + self.sC_0 - self.sC_3)) - (
                         2 * self.x0 + A[0] + A[1] + A[2] + A[3] + A[5]),
                 int(0.5 * (
                             - self.sD_0 + self.sD_1 - 2 * self.sR_0 - self.sR_1 - self.sR_2 + self.sC_0 + self.sC_3)) + (
                         self.x0 + A[0] + A[1] + A[3] + A[5]),
                 int(0.5 * (
                             - self.sD_0 - self.sD_1 - self.sR_1 - self.sR_2 + 2 * self.sC_1 + self.sC_0 + self.sC_3)) + (
                         -A[0] + A[4] + A[5]),
                 self.sC_2 - (A[1] + A[4] + A[5]))
            # Since all elements of the magic squares are all different, then could be: B = (A u {x0})^c = C
            if set(B) == C:
                # Build a square with ordered element:
                new_square = (self.x0, A[0], A[1], B[0], A[2], A[3], A[4], B[1], B[5], B[4], A[5], B[3], B[6],
                             B[7], B[8], B[2])
                if show:
                    print(new_square)
                squares.add(new_square)
        # Delete repeated isomorphic squares
        squares = reduce_isomorphic_squares(squares)
        # Save squares:
        time_now = datetime.now().strftime("%y.%m.%d_%H.%M.%S")
        if save:
            with open(output_directory + "/squares_type_I_" + time_now + ".txt", "w") as f:
                k = 0
                format_k = '{:>' + str(int(np.log10(len(squares)+1.00000001)) + 1) + '}'
                for M in squares:
                    k += 1
                    # With label
                    # f.write(format_k.format(k) + ' ' + str(M) + '\n')
                    # Without label
                    f.write(str(M) + '\n')
        return squares

    def compute_squares_of_type_II(self, show=False, save=False):
        # Magic square of type I
        # A[0]  x0  A[1] B[0] \\
        # A[2] A[3] A[4] B[1] \\
        # B[5] B[4] A[5] B[3] \\
        # B[6] B[7] B[8] B[2]
        if show:
            print("Computing squares of Type II \n", flush=True)
        # Remove x0 from the set of element possibilities
        domain = self.elements.difference({self.x0})
        # Compute all combinations of 6 elements taken from domain (15 elements)
        comb = permutations(domain, 6)
        # Output set. This set contains also no ordinary magic squares.
        squares = set()
        for A in comb:
            # Compute the complement of A u {x0} (9 elements)
            C = domain - set(A)
            # Compute the 9 dependent variables
            B = (self.sR_0 - (self.x0 + A[0] + A[1]),
                 self.sR_1 - (A[2] + A[3] + A[4]),
                 self.sD_0 - (A[0] + A[3] + A[5]),
                 self.sC_3 - (self.sR_0 + self.sR_1 + self.sD_0) + (
                             self.x0 + 2*A[0] + A[1] + A[2] + 2 * A[3] + A[4] + A[5]),
                 int(0.5 * (self.sD_0 + self.sD_1 + self.sR_1 + self.sR_2 - self.sC_0 - self.sC_3)) - (
                             A[3] + A[4] + A[5]),
                 int(0.5 * (self.sD_0 - self.sD_1 + 2 * self.sR_0 + self.sR_1 + self.sR_2 + self.sC_0 - self.sC_3)) - (
                         self.x0 + 2*A[0] + A[1] + A[2] + A[3] + A[5]),
                 int(0.5 * (
                             - self.sD_0 + self.sD_1 - 2 * self.sR_0 - self.sR_1 - self.sR_2 + self.sC_0 + self.sC_3)) + (
                         self.x0 + A[0] + A[1] + A[3] + A[5]),
                 int(0.5 * (
                             - self.sD_0 - self.sD_1 - self.sR_1 - self.sR_2 + 2 * self.sC_1 + self.sC_0 + self.sC_3)) + (
                         -self.x0 + A[4] + A[5]),
                 self.sC_2 - (A[1] + A[4] + A[5]))
            # Since all elements of the magic squares are all different, then could be: B = (A u {x0})^c = C
            if set(B) == C:
                # Build a square with ordered element:
                new_square = (A[0], self.x0, A[1], B[0], A[2], A[3], A[4], B[1], B[5], B[4], A[5], B[3], B[6],
                             B[7], B[8], B[2])
                if show:
                    print(new_square)
                squares.add(new_square)
        # Delete repeated isomorphic squares
        squares = reduce_isomorphic_squares(squares)
        # Save sorted squares:
        time_now = datetime.now().strftime("%y.%m.%d_%H.%M.%S")
        if save:
            with open(output_directory + "/squares_type_II_" + time_now + ".txt", "w") as f:
                k = 0
                format_k = '{:>' + str(int(np.log10(len(squares)+1.00000001)) + 1) + '}'
                for M in squares:
                    k += 1
                    # With label
                    # f.write(format_k.format(k) + ' ' + str(M) + '\n')
                    # Without label
                    f.write(str(M) + '\n')
        return squares

    def compute_squares_of_type_III(self, squares_type_I, show=False, save=False):
        if show:
            print("Computing squares of Type III \n", flush=True)
        squares = set()
        for M in squares_type_I:
            new_square = permutation_I_to_III(M)
            if show:
                print(new_square)
            squares.add(new_square)
        # Save sorted squares:
        time_now = datetime.now().strftime("%y.%m.%d_%H.%M.%S")
        if save:
            with open(output_directory + "/squares_type_III_" + time_now + ".txt", "w") as f:
                k = 0
                format_k = '{:>' + str(int(np.log10(len(squares)+1.00000001)) + 1) + '}'
                for M in squares:
                    k += 1
                    # With label
                    # f.write(format_k.format(k) + ' ' + str(M) + '\n')
                    # Without label
                    f.write(str(M) + '\n')
        return squares

    def compute_magic_squares(self, show=False, save=False):
        squares = self.compute_squares_of_type_I(show=show,save=save)
        squares = squares.union(self.compute_squares_of_type_III(squares,show=show,save=save))
        squares = squares.union(self.compute_squares_of_type_II(show=show,save=save))
        print(f"TOTAL OF CLASES MAGIC SQUARES = {len(squares)}", flush=True)

    def compute_hiper_magic_squares(self, show = False, save = False):
        squares = self.compute_squares_of_type_I()
        squares = squares.union(self.compute_squares_of_type_III(squares))
        squares = squares.union(self.compute_squares_of_type_II())
        hiper_magic_squares = set()
        for M in squares:
            if is_hiper_magic(M,self.s):
                if show:
                    print(M)
                hiper_magic_squares.add(M)
        # Save squares:
        time_now = datetime.now().strftime("%y.%m.%d_%H.%M.%S")
        if save:
            with open(output_directory + "/hiper_magic_squares_" + time_now + ".txt", "w") as f:
                k = 0
                format_k = '{:>' + str(int(np.log10(len(hiper_magic_squares))) + 1) + '}'
                for M in hiper_magic_squares:
                    k += 1
                    # With label
                    # f.write(format_k.format(k) + ' ' + str(M) + '\n')
                    # Without label
                    f.write(str(M) + '\n')
        print(f"TOTAL OF CLASES OF HIPER MAGIC SQUARES = {len(hiper_magic_squares)}", flush=True)

    def compute_super_magic_squares(self, show = False, save = False):
        squares = self.compute_squares_of_type_I()
        squares = squares.union(self.compute_squares_of_type_III(squares))
        squares = squares.union(self.compute_squares_of_type_II())
        super_magic_squares = set()
        for M in squares:
            if is_super_magic(M,self.s):
                if show:
                    print(M)
                super_magic_squares.add(M)
        # Save squares:
        time_now = datetime.now().strftime("%y.%m.%d_%H.%M.%S")
        if save:
            with open(output_directory + "/super_magic_squares_" + time_now + ".txt", "w") as f:
                k = 0
                format_k = '{:>' + str(int(np.log10(len(super_magic_squares))) + 1) + '}'
                for M in super_magic_squares:
                    k += 1
                    # With label
                    # f.write(format_k.format(k) + ' ' + str(M) + '\n')
                    # Without label
                    f.write(str(M) + '\n')
        print(f"TOTAL OF CLASES OF SUPER MAGIC SQUARES = {len(super_magic_squares)}", flush=True)

def reduce_isomorphic_squares(squares):
    # Delete repeated isomorphic squares under rotations and reflections
    new_squares = set()
    for M in squares:
        k = 0
        if rotate_90(M) not in new_squares:
            k+=1
        if rotate_180(M) not in new_squares:
            k+=1
        if rotate_270(M) not in new_squares:
            k += 1
        if reflect_0(M) not in new_squares:
            k += 1
        if reflect_45(M) not in new_squares:
            k += 1
        if reflect_90(M) not in new_squares:
            k += 1
        if reflect_135(M) not in new_squares:
            k += 1
        if k==7:
            new_squares.add(M)
    return new_squares

def is_hiper_magic(M,s):
    # M[0]  M[1]  M[2]  M[3] \\
    # M[4]  M[5]  M[6]  M[7] \\
    # M[8]  M[9]  M[10] M[11] \\
    # M[12] M[13] M[14] M[15]
    if M[1] + M[4] + M[11] + M[14] == s and M[2] + M[7] + M[8] + M[13] == s :
        return True
    else:
        return False

def is_super_magic(M,s):
    # M[0]  M[1]  M[2]  M[3] \\
    # M[4]  M[5]  M[6]  M[7] \\
    # M[8]  M[9]  M[10] M[11] \\
    # M[12] M[13] M[14] M[15]
    if M[4] + M[5] + M[8] + M[9] == s and M[6] + M[7] + M[10] + M[11] == s and \
            M[1] + M[2] + M[5] + M[6] == s and M[9] + M[10] + M[13] + M[14] == s:
        return True
    else:
        return False

def rotate_90(M):
    # M is a vector of 16 elements
    # Rotate 90 degrees anticlockwise
    return (M[3],M[7],M[11],M[15],M[2],M[6],M[10],M[14],M[1],M[5],M[9],M[13],M[0],M[4],M[8],M[12])

def rotate_180(M):
    # M is a vector of 16 elements
    # Rotate 180 degrees anticlockwise
    return (M[15],M[14],M[13],M[12],M[11],M[10],M[9],M[8],M[7],M[6],M[5],M[4],M[3],M[2],M[1],M[0])

def rotate_270(M):
    # M is a vector of 16 elements
    # Rotate 270 degrees anticlockwise
    return (M[12],M[8],M[4],M[0],M[13],M[9],M[5],M[1],M[14],M[10],M[6],M[2],M[15],M[11],M[7],M[3])

def reflect_90(M):
    # M is a vector of 16 elements
    # Reflect respect Y axis
    return (M[3],M[2],M[1],M[0],M[7],M[6],M[5],M[4],M[11],M[10],M[9],M[8],M[15],M[14],M[13],M[12])

def reflect_0(M):
    # M is a vector of 16 elements
    # Reflect respect X axis or Rot(180)Ref(90)
    return (M[12],M[13],M[14],M[15],M[8],M[9],M[10],M[11],M[4],M[5],M[6],M[7],M[0],M[1],M[2],M[3])

def reflect_135(M):
    # M is a vector of 16 elements
    # Reflect respect Y=-X or Rot(90)Ref(90)
    return (M[0],M[4],M[8],M[12],M[1],M[5],M[9],M[13],M[2],M[6],M[10],M[14],M[3],M[7],M[11],M[15])

def reflect_45(M):
    # M is a vector of 16 elements
    # Reflect respect Y=X or Rot(270)Ref(90)
    return (M[15],M[11],M[7],M[3],M[14],M[10],M[6],M[2],M[13],M[9],M[5],M[1],M[12],M[8],M[4],M[0])

def transposition_row_02(M):
    # M is a vector of 16 elements
    return (M[8],M[9],M[10],M[11],M[4],M[5],M[6],M[7],M[0],M[1],M[2],M[3],M[12],M[13],M[14],M[15])

def transposition_col_02(M):
    # M is a vector of 16 elements
    return (M[2],M[1],M[0],M[3],M[6],M[5],M[4],M[7],M[10],M[9],M[8],M[11],M[14],M[13],M[12],M[15])

def transposition_row_13(M):
    # M is a vector of 16 elements
    return (M[0],M[1],M[2],M[3],M[12],M[13],M[14],M[15],M[8],M[9],M[10],M[11],M[4],M[5],M[6],M[7])

def transposition_col_13(M):
    # M is a vector of 16 elements
    return (M[0],M[3],M[2],M[1],M[4],M[7],M[6],M[5],M[8],M[11],M[10],M[9],M[12],M[15],M[14],M[13])

def diagonal_permutation_02(M):
    return transposition_col_02(transposition_row_02(M))

def diagonal_permutation_13(M):
    return transposition_col_13(transposition_row_13(M))

def permutation_I_to_III(M):
    return reflect_0(diagonal_permutation_02(diagonal_permutation_13(reflect_90(M))))


