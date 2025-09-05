#!/home/joonho/anaconda3/bin/python

import argparse
import re
import os
from common import f_ext
import openpyxl

def extract_excel(xlsx, ext):
    # Load the workbook
    workbook = openpyxl.load_workbook(xlsx)

    # Select the active sheet (or use workbook['SheetName'])
    sheet = workbook.active

    # Example: read all rows
    for row in sheet.iter_rows(values_only=True):
        print(row)
    return 0

def main():
    parser = argparse.ArgumentParser(description='To get rid of abnormal DOS at start energy')
    parser.add_argument('file', nargs='?', default='DOSCAR', help='read DOSCAR')
    parser.add_argument('-o', '--old', default='_o', help='save original DOSCAR to DOSCAR_o')
    args = parser.parse_args()

    if f_ext(args.file) != 'xlsx':
        print(f"Warning: file is not extension of xlsx")
    
    extract_excel(args.file, args.old) 

if __name__ == '__main__':
    main()
