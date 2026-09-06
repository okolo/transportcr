#!/usr/bin/python3
'''
Created on Jan 16, 2014

@author: ok
'''

import sys

from xsw_utils import (
    XswParameterError,
    replace_parameter_in_text,
    validate_parameter_value,
    write_text_atomic,
)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv) < 4:
        print('usage: replacePar <parName> <parValue> <file1.xsw> [file2.xsw [..]]')
        return 1

    try:
        newValue = validate_parameter_value(argv[2])
    except XswParameterError as error:
        print('error: ' + str(error), file=sys.stderr)
        return 1
    for iArg in range(3, len(argv)):
        xswFileName = argv[iArg]
        with open(xswFileName, "r") as xswFile:
            original = xswFile.read()
        updated, nReplacements = replace_parameter_in_text(
            original, argv[1], newValue
        )
        if nReplacements == 0:
            print(xswFileName + ': parameter not found')
        else:
            if nReplacements > 1:
                tmpFileName = xswFileName + ".tmp"
                with open(tmpFileName, "w") as tmpFile:
                    tmpFile.write(updated)
                print(xswFileName + ': ' + str(nReplacements) + ' replacements were made!!! Original file unchanged. Modified file saved in ' + tmpFileName)
            else:
                write_text_atomic(xswFileName, updated)
    return 0


if __name__ == '__main__':
    sys.exit(main())
