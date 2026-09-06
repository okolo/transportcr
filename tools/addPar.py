#!/usr/bin/env python3
'''
Created on Jan 16, 2014

@author: ok
'''

import sys
import os


def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv) < 4:
        print('usage: addPar <parName> <sampleFile.xsw> <file1.xsw> [file2.xsw [..]]')
        print('if param already exists in file1.xsw etc. it will be duplicated! No checks are performed!')
        return 1

    paramName=argv[1]
    titlePattern = "title=\"" + paramName  + "\""
    anyTitlePattern = "title=\""
    startTagPattern = "<"
    endTagPattern = "/>"
    valuePattern = "value=\""
    sampleFileName = argv[2]
    tagLine=""
    tagFound = False
    nextParName=""
    with open(sampleFileName, "r") as sampleFile:
        for line in sampleFile.readlines():
            if tagFound:
                j=line.find(valuePattern)
                if j<=0:
                    continue
                i = line.find(anyTitlePattern)
                nextParName = line[i + len(anyTitlePattern):line.find("\"", i + len(anyTitlePattern))]
                break
            i=line.find(titlePattern)
            if i>=0:
                iTagStart=line.rfind(startTagPattern, 0, i)
                iTagEnd=line.find(endTagPattern, i)
                if iTagStart>=0 and iTagEnd>0:
                    tagLine = line[iTagStart:iTagEnd+len(endTagPattern)]
                    tagFound=True

    if not tagFound:
        print("parameter " + paramName + " not found in " + sampleFileName)
        return 1
    if len(nextParName)==0:
        nextParName="MissingParamsCount"

    nextParPattern = "title=\"" + nextParName  + "\""
    for iArg in range(3, len(argv)):
        xswFileName = argv[iArg]
        with open(xswFileName, "r") as xswFile:
            tmpFileName = xswFileName + ".tmp"
            with open(tmpFileName, "w") as xswFile2:
                added=False
                for line in xswFile.readlines():
                    if len(line.strip()) == 0:  # print only nonempty lines
                        continue
                    if added:
                        xswFile2.write(line)
                        continue
                    i=line.find(nextParPattern)
                    if i>=0:
                        iNextTagStart=line.rfind(startTagPattern, 0, i)
                        if iNextTagStart>=0:
                            xswFile2.write(line[0:iNextTagStart] + tagLine + "\n")
                            nSpaces = len(line)-len(line.lstrip())
                            xswFile2.write(line[0:nSpaces] + line[iNextTagStart:])
                            added=True
                        else:
                            print("multiline neighbour tags not supported, " + xswFileName + " unchanged, see " + tmpFileName)
                            return 1
                    else:
                        xswFile2.write(line)

                if not added:
                    print("neighbour parameter " + nextParName + ' not found in ' + xswFileName + ", file unchanged")
                    os.remove(tmpFileName)
                else:
                    os.remove(xswFileName)
                    os.rename(tmpFileName, xswFileName)
    return 0


if __name__ == '__main__':
    sys.exit(main())
