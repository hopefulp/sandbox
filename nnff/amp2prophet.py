#!/home/joonho/anaconda3/bin/python

import argparse
from amp import Amp
from amp.convert import save_to_prophet

def make_prophet(amppot):
    calc = Amp.load(amppot)
    save_to_prophet(calc)

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--amppot', default='amp.amp', help='amp potential')
    args = parser.parse_args()

    make_prophet(args.amppot)

