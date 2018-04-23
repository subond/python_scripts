from netcdftime import utime
from datetime import  datetime
from cmip_time import FakeDT
import numpy as np

def day_number_to_datetime_array(time_in, calendar_type, units_in):

	cdftime = utime(units_in, calendar = calendar_type)

	date_out = cdftime.num2date(time_in)

	return date_out

def day_number_to_date(time_in, calendar_type = '360_day', units_in = 'days since 0001-01-01 00:00:00'):
	"""
	Aim is to make the time array have attributes like .month, or .year etc. This doesn't work with
	normal datetime objects, so Mike's FakeDT does this for you. First step is to turn input times
	into an array of datetime objects, and then FakeDT makes the array have the attributes of the
	elements themselves.
	"""

	time_in = day_number_to_datetime_array(time_in, calendar_type, units_in)


	cdftime = FakeDT( time_in, units=units_in,
                 calendar=calendar_type)

	return cdftime

def month_to_season(months_in):

	seasons=np.zeros(len(months_in))
	idx_djf=(months_in==1)&(months_in==2)&(months_in==12)	
	idx_mam=(months_in==3)&(months_in==4)&(months_in==5)
	idx_jja=(months_in==6)&(months_in==7)&(months_in==8)
	idx_son=(months_in==9)&(months_in==10)&(months_in==11)

	seasons[idx_djf]=0
	seasons[idx_mam]=1
	seasons[idx_jja]=2
	seasons[idx_son]=3


	return seasons




if __name__ == "__main__":
    import numpy as np
    from datetime import  datetime

    cdftime = utime('hours since 0001-01-01 00:00:00')
    date = datetime.now()
    print date
    t = cdftime.date2num(date)
    print t
    date = cdftime.num2date(t)
    print date

    time = 1772.5
    date_equivalent = day_number_to_date(time)
    print date_equivalent





