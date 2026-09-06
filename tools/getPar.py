#!/usr/bin/python3
'''
Created on Jan 16, 2014

@author: ok
'''

import sys

from xsw_utils import get_parameter_values

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv) != 3:
        print('usage: getPar <parName> <file.xsw>')
        return 1

    for value in get_parameter_values(argv[2], argv[1]):
        print(value)
    return 0


if __name__ == '__main__':
    sys.exit(main())
