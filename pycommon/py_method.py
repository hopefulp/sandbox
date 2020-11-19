'''
    special methods
'''

import re
import subprocess


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def get_cli(com):
    '''
    return shell command output
    '''
    #com = "grep 'optimization ' %s/amp-log.txt -B 1 | awk '/T/ {print $8, $10}'" % x
    results = subprocess.check_output(com, shell=True)
    #yv = float(results.split()[0])
    #y.append(yv)
    return results

