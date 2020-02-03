from __future__ import print_function
import lemon
import os

file1 = lemon.open_file("../../../test/files/1AAQ.mmtf")
file2 = lemon.open_file("../../../test/files/1YT9.mmtf.gz")

transformation = lemon.kabsch(file1.positions_const(), file2.positions_const(), 1e-10)

my_pos = file2.positions()

lemon.align(my_pos, transformation)

lemon.write_file(file2, "out.pdb")
os.remove("out.pdb")
