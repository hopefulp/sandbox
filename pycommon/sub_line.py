#!/home/joonho/anaconda3/bin/python

import server_env


def sub_line_dict(fname, dic):
    with open(fname, 'rw') as f:
        lines = f.readlines()
        for key in dic.keys():
            sub_line(lines, key, dic['key'])
    return 0            


def sub_line_dict(fname, jobtype):
    return 0

def main():
    parser = argparse.ArgumentParser(description='substitution a line')
    parser.add_argument('fname', help='input file')
    parser.add_argument('jobtype', choices=['python'], help='job type')
    args = parser.parse_args()

    sub_line(args.fname, args.jobtype) 
    return 0

if __name__ == '__main__':
    main()
