#!/home/joonho/anaconda3/bin/python

import argparse

d_search=['/home/joonho/bin', '/home/joonho/sandbox_gl']

def jobs(Lclass, job, Lusage):


    return 0

def main():
    parser = argparse.ArgumentParser(description="display Usage for /mymplot  ")
    #parser.add_argument('-c','--classify', action='store_true', help="classify ")
    #parser.add_argument('-j','--job', help="present class details ")
    #parser.add_argument('-js','--specify', choices=['qcmo','nbo','eda'], help="present class details ")
    #parser.add_argument('-u','--usage', action='store_true', help="present main details")
    args = parser.parse_args()

    jobs(args.classify, args.job, args.usage)

if __name__ == "__main__":
    main()
