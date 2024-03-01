import re
import statistics as st


MYREGEX = " time cost = *(\d*\.[0-9]*)"

log_path = f'./run_logs/wbctry.log'

log_file = open(log_path, "r")
log_info = log_file.read()
log_file.close()

times_list = list(map(float, re.findall(MYREGEX, log_info)))

median_time = st.median(times_list)

ostr = f"Median time: {median_time}"
print(ostr)