#!/bin/env python

from os import popen

#numero iniziale
i=197273588
  

#inizio ciclo
while i<197274095:
    out = popen(" echo \"bkill "+str(i)+" \" ")
    for x in out:
        print x
    out=popen("bkill "+str(i)+"" )
    for x in out:
        print x
    i+=1
