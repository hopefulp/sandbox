import os

os.system('nebmake.pl POSCAR POSCAR 1')
os.system('rm -rf 01 02')
os.system('mv 00/POSCAR .')
os.system('rm -rf 00')
