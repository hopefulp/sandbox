#!/home/noische/program/python27/bin/python 
"""
smartMerge.py
Original: Jan 25 2011 In Kim

find a cross-section of two text files and merge them

usage: smartMerge.py filename1 filename2 output
"""

import sys
import string

version = '110125'

def long_substr(data):
    """
long_substr(data):
    find the longest common substring between the elements of data
    """
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and is_substr(data[0][i:i+j], data):
                    substr = data[0][i:i+j]
    return substr


def is_substr(find, data):
    """
is_substr(find, data):
    query whether a string is a subset of a data or not
    """
    if len(data) < 1 and len(find) < 1:
        return False
    for i in range(len(data)):
        if find not in data[i]:
            return False
    return True


def readFileInLine(filename):
    """
readFileInLine(filename):
    reads a file and returns one string
    whole file contents are in a line
    """
    try:
        filehandler = open(filename)
        lines = filehandler.readlines()
    except IOError:
        print("Error: Filename " + filename + " cannot be read.")
    else:
        filehandler.close()

    line = "".join(lines)
    line = line.replace("\r\n", "")

    return line


def reverseString(str):
    """
reverseString(str):
    returns a string which is reversed by characters
    """
    revstr = list(str)
    revstr.reverse()
    revstr = "".join(revstr)

    return revstr


if __name__ == '__main__':

    a = readFileInLine(sys.argv[1])	# read the first file
    b = readFileInLine(sys.argv[2])	# read the second file

    result_file = open(sys.argv[3], 'w')	# prepare output file
    error = ""
    
    a = a.replace("N","")	# remove all N's in the file
    b = b.replace("N","")	# remove all N's in the file

    cross = long_substr((a, b))	# find the cross-section of two file
    a1, a2, a3 = a.partition(cross)	# the first file is split into pure a, cross-section, and NULL. notice that a2 == b2.
    b1, b2, b3 = b.partition(cross)	# the second file is split into NULL, cross-section, and pure b

    if len(b1) != 0:	# if the length of the cross-section is too short, print warning.
        error = "Warning: " + sys.argv[1] + " and " + sys.argv[2] + " seems not to be fully matched. You may check the input files.\n"
        print(error)

    result = []			# temporary variable
    result.append(a1)		# a contains pure 'a' part and cross-section
    result.append(" " + cross + " ")
    result.append(b3)		# result = ('pure a', 'cross-section', 'pure b')
    result = "".join(result)	# concatenate result into string
    result += "\n"

    output = ""
    output += "*** Contents in " + sys.argv[1] + "\n"
    output += a1 + " " + a2 + " " + a3 + "\n\n"
    output += "*** Contents in " + sys.argv[2] + "\n"
    output += b1 + " " + b2 + " " + b3 + "\n\n"
    output += "*** Common substring: \n"
    output += cross + "\n\n"
    output += error
    output += "*** Result: \n"
    output += result + "\n"

    result_file.write(output)	# write to the file
