#!/usr/bin/python3
'''
Created on Jan 16, 2014

@author: ok
'''

import sys
from pathlib import Path


def default_xsw_path():
    return Path(__file__).resolve().parent.parent / "bin" / "switches.xsw"


def main(argv=None):
    if argv is None:
        argv = sys.argv
    if len(argv) < 2 or argv[1] == "--help":
        print('usage: ' + argv[0] + ' <parName> [file.xsw]')
        print("default xsw file is " + str(default_xsw_path()))
        return 1

    titlePatern = "title=\"" + argv[1] + "\""
    commentPatern= "comment=\""
    xswFileName = str(default_xsw_path())
    if len(argv) > 2:
        xswFileName = argv[2]
    xswFile = open(xswFileName, "r")
    switchPrint = False
    curSw=0
    for line in xswFile.readlines():
        line=line.strip()
        totLen=len(line)
        if switchPrint:
            iEndSw=line.find("</switch>")
            if iEndSw>=0:
                totLen=iEndSw
            iStart=0
            while iStart<totLen:
                iNextOption=line.find("<sw",iStart)
                if iNextOption<0:
                    break
                iEndOption=line.find("/>",iNextOption)
                if iEndOption<0:
                    iEndOption=totLen-2
                print(str(curSw) + "\t" + line[iStart:iEndOption+2])
                curSw += 1
                iStart=iEndOption+2
            if iEndSw>=0:
                break
        i=line.find(titlePatern)
        if i>=0:
            iCommentStart=line.find(commentPatern)
            if iCommentStart>=0:
                iCommentStart = iCommentStart + len(commentPatern)
                iValEnd=line.find("\"", iCommentStart)
                if iValEnd>iCommentStart:
                    print(line[iCommentStart:iValEnd])
            if line.rfind("<switch",0,i)>=0: #if par is switch type
                switchPrint = True

    xswFile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
