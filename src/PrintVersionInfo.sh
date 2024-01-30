#!/bin/bash
ver=`git log --oneline | wc -l | awk ' { printf("%d", $1-866+1000) }'`
modified=`git diff . | wc -l`
if [ $modified -gt 0 ]
then
	ver=${ver}M
fi
hash=`git show -s --pretty=format:%H`
echo -n "$ver : $hash"

