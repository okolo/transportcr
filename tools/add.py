#!/usr/bin/env python3
'''
Created on Dec 3, 2010

@author: kalashev
'''
from __future__ import print_function
import sys
import math

replace_negative_by_zero=True

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def tabfunc(x, table, yCol, xCol=0):
    xData=table[xCol]

    if x<xData[0]:
        return 0
    i2=len(xData)-1
    if x>xData[i2]:
        return 0
    i1=0

    while i2-i1>1:
        i=int((i1+i2)/2)
        if(x>xData[i]):
            i1=i
        else:
            i2=i
    x1=xData[i1]
    x2=xData[i2]
    y1=table[yCol][i1]
    y2=table[yCol][i2]
    return y1 + (x-x1)*(y2-y1)/(x2-x1)

def readTable(fileName):
    result = []
    lineNo = 0
    with open(fileName, 'r') as f:
        for line in f:
            lineNo = lineNo+1
            stripped_line = line.strip()
            if not stripped_line or stripped_line.startswith('#'):
                continue
            dline = stripped_line.split()
            if len(result) == 0:
                for i in range(len(dline)):
                    result.append([float(dline[i])])
            else:
                if len(dline) != len(result):
                    eprint('input format error: variable number of columns in ' + fileName + ' line ' + str(lineNo))
                    sys.exit(1)
                for i in range(len(dline)):
                    result[i].append(float(dline[i]))
    return result

def buildScale(tables, xCol=0):
    result = []
    minX=1e300
    maxX=-1e-300
    minStep=1e300
    for data in tables:
        xMin = data[xCol][0]
        xMax = data[xCol][len(data[xCol])-1]
        step = math.log10(xMax/xMin)/(len(data[xCol])-1)
        if xMin<minX:
            minX=xMin
        if xMax>maxX:
            maxX=xMax
        if step<minStep:
            minStep=step
    step=math.pow(10.0, step)
    x=minX
    while x<=maxX:
        result.append(x)
        x=x*step
    return result

def add(multipliers, components, output=sys.stdout, zero_negative=None):
    assert len(multipliers) == len(components)
    if zero_negative is None:
        zero_negative = replace_negative_by_zero
    data = [readTable(f) for f in components]

    nCols = len(data[0])

    for nComp in range(0, len(data)):
        if len(data[nComp]) != nCols:
            eprint('component ' + str(nComp + 1) + ' has ' + str(
                len(data[nComp])) + ' columns while component 0 has ' + str(nCols) + ' columns')
            sys.exit(1)

    for x in buildScale(data):
        line = str(x)
        for col in range(1, nCols):
            val = 0.0
            for nComp in range(0, len(data)):
                val = val + (multipliers[nComp]) * tabfunc(x, data[nComp], col)
            if zero_negative and val < 0:
                val = 0.
            line = line + '\t' + '{0:.16g}'.format(val)
        print(line, file=output)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv)<3 or len(argv)%2==0:
        eprint('Usage: add.py multiplier1 component1 [multiplier2 component2 [...]] ')
        return 1

    files = []
    mult = []

    for nf in range(2,len(argv),2):
        mult.append(float(argv[nf-1]))
        files.append(argv[nf])

    add(mult, files)
    return 0


if __name__ == '__main__':
    sys.exit(main())
