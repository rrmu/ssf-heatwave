# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 21:05:52 2022

@author: Michael

EXECUTION PROGRAM FOR HEATWAVE ANALYSIS
"""
import tempcalc as tc
import exporter as exp
import os
import concurrent.futures

reffile = 'D:/Temperature/tmax_climatology_8110.nc'

print("===PREPARING DATA===")

####GET DIRECTORY LIST####
def completedir(dirname):
    newdir=[]
    for file in os.listdir(dirname):
        newdir.append(os.path.join(dirname,file).replace("\\","/"))
    return newdir
        
dir1=completedir('D:/Temperature/1')
dir2=completedir('D:/Temperature/2')
dir3=completedir('D:/Temperature/3')
dir4=completedir('D:/Temperature/4')

####GET PERCENTILE####

percentile=tc.tc.percentile(reffile, .9)
"""
print("===ANALYZING HEATWAVE CASE===")

####CALCULATE YEARLY HEAT WAVE OCCURENCE####
yearbool8089=tc.multithread.daysTF(dir1,p1,3)
yearbool9099=tc.multithread.daysTF(dir2,p2,3)
yearbool0009=tc.multithread.daysTF(dir3,p3,3)
yearbool1022=tc.multithread.daysTF(dir4,p4,3)
print()
print()
print('!!Bool List Ready!!')
print()
print()


####CALCULATE YEARLY HEAT WAVE EVENT COUNT####
eventpile=[]

with concurrent.futures.ProcessPoolExecutor(4) as multiprocess:
    event_futures=[multiprocess.submit(tc.multithread.daysTF,DIR,percentile,3,isevent=True)
                   for DIR in [dir1,dir2,dir3,dir4]]

    eventpile.append(ten_year_futures.result() 
                     for ten_year_futures in concurrent.futures.as_completed(event_futures))

yearevent=[]
eventfin=[]
for yearlists in eventpile:
    eventfin.append(yearlists.result())
    for yeartup in eventfin:
        yearevent.append(yeartup)

yearevent=sorted(yearevent)

print()
print()
print('!!Event List Ready!!')
print()
print()


"""

yearevent8089=tc.multithread.daysTF(dir1,percentile,3,isevent=True)
yearevent9099=tc.multithread.daysTF(dir2,percentile,3,isevent=True)
yearevent0009=tc.multithread.daysTF(dir3,percentile,3,isevent=True)
yearevent1022=tc.multithread.daysTF(dir4,percentile,3,isevent=True)




"""

####CALCULATE YEARLY HEAT WAVE DAY COUNT####
year8089=tc.multithread.puredaycount(yearbool8089)
year9099=tc.multithread.puredaycount(yearbool9099)
year0009=tc.multithread.puredaycount(yearbool0009)
year1022=tc.multithread.puredaycount(yearbool1022)

print()
print()
print('!!Day List Ready!!')
print()
print()

"""



yearevent=[]
for tup in yearevent8089:
    yearevent.append(tup)
    
for tup in yearevent9099:
    yearevent.append(tup)
    
for tup in yearevent0009:
    yearevent.append(tup)
   
for tup in yearevent1022:
    yearevent.append(tup)


####DATA PLOTTING####




"""
####EXPORT PERCENTILES####
exp.pt2c(reffile1,p1,'percentile1')
exp.pt2c(reffile2,p2,'percentile2')
exp.pt2c(reffile3,p3,'percentile3')
exp.pt2c(reffile4,p4,'percentile4')

####EXPORT DATA####
exp.cdf.event10Yr(reffile1,yearevent8089,'hwcase','HeatWaveOccurence')
exp.cdf.event10Yr(reffile2,yearevent9099,'hwcase','HeatWaveOccurence')
exp.cdf.event10Yr(reffile3,yearevent0009,'hwcase','HeatWaveOccurence')
exp.cdf.event10Yr(reffile4,yearevent1022,'hwcase','HeatWaveOccurence')
"""