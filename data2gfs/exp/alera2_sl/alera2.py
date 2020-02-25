#!/usr/bin/env python
from datetime import datetime
from dateutil.relativedelta import relativedelta
import os
import os.path
from subprocess import call
import sys

date=datetime(2009,9,30,12)
os.environ['JCAP']='190'
os.environ['IDATE']=date.strftime('%Y%m%d%H')

for member in range(64):
    if member==0:
        os.environ['MEMBER']='mean'
    else:
        os.environ['MEMBER']='%03d'%member
    ret = call(['./alera2_sl.sh'])
    if ret!= 0:
        print ret
        break

    
