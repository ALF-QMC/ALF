#!/bin/sh
if git status  > /dev/null 2>&1
then 
  echo "#define GIT" > git.h
  if [ -n "$(git status --porcelain > /dev/null 2>&1)" ]; then
    echo "#define GIT_COMMIT_HASH \"$(git log -1 --format=%h)-dirty\"" >> git.h
  else
    echo "#define GIT_COMMIT_HASH \"$(git log -1 --format=%h)\"" >> git.h
  fi
  echo "#define GIT_BRANCH \"$(git rev-parse --abbrev-ref HEAD)\"" >> git.h
  git status | sed 's|"|""|g' | sed 's|^|write(50,*) "|' | sed 's|$|"|' > git_status.h
else
  echo "" > git.h
  echo "" > git_status.h
fi
