import datetime as dt
import os

for yy in range(1980,2020):
	for mm in range(1,13):
		#print(yy, "{:02d}".format(mm))
		date = dt.date(yy, mm, 1)
		mm_ = mm
		while mm_ == mm:
			print(date)
			day_ = str(date.day)
			mm_ = str(mm)
			if len(day_)==1:
				day_ = "0"+day_
			if len(mm_)==1:
				mm_ = "0"+mm_
			print("wget https://cfs.ncep.noaa.gov/cfs/godas/pentad/"+str(yy)+"/godas.P."+str(yy)+mm_+day_+".grb")
			cmd = "wget https://cfs.ncep.noaa.gov/cfs/godas/pentad/"+str(yy)+"/godas.P."+str(yy)+mm_+day_+".grb"
			os.system(cmd)
			date = date + dt.timedelta(days=1)
			mm_ = date.month
