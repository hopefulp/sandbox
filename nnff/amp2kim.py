#!/home/joonho/anaconda3/bin/python

import argparse


def jobs(Lclass, job, Lusage):


    return 0

def main():
    parser = argparse.ArgumentParser(description="  ")
    #parser.add_argument('-c','--classify', action='store_true', help="classify ")
    #parser.add_argument('-u','--usage', action='store_true', help="present main details")
    args = parser.parse_args()

    jobs(args.classify, args.job, args.usage)

if __name__ == "__main__":
    main()
