#!/bin/sh
if [ -d ../.git ]
then 
  echo "#define GIT" > git.h
  echo "#define GIT_COMMIT_HASH \"$(git log -1 --format=%h)\"" >> git.h
  echo "#define GIT_BRANCH \"$(git rev-parse --abbrev-ref HEAD)\"" >> git.h
else
  echo "" > git.h
fi
