import re
import statistics as st


MYREGEX = " time cost = *(\d*\.[0-9]*)"

print ("Procs, Layers, Time")

for i in range(1,6):
    l_ct = 2**i
    proc_ct = l_ct * 6
    log_path = f'./run_logs/250I_{l_ct}L_{proc_ct}P.log'

    log_file = open(log_path, "r")
    log_info = log_file.read()
    log_file.close()

    times_list = list(map(float, re.findall(MYREGEX, log_info)))

    median_time = st.median(times_list)

    ostr = f'{proc_ct}, {l_ct}, {median_time}'
    print(ostr)