"""
from http://stackoverflow.com/questions/5704931/parse-string-of-integer-sets-with-intervals-to-list
>>> s = "2,5,7-9,12"
>>> ranges = (x.split("-") for x in s.split(","))
>>> print [i for r in ranges for i in range(int(r[0]), int(r[-1]) + 1)]
[2, 5, 7, 8, 9, 12]
"""

