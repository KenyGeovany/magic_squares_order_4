from fast_general_magic_squares_4x4 import *

M = magic_square()
M.elements = {1,2,3,4,5,6,7,8,9,11,13,15,16,18,20,22}
M.s = 40
M.set_isotropic()
M.compute_magic_squares(show=True, save=True)