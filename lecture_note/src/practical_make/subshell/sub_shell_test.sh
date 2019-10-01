#!/bin/sh

echo My process id and shell level = $$ $SHLVL

if [ $SHLVL -gt 3  ]
then
  exit
fi

sh ./sub_shell_test.sh

