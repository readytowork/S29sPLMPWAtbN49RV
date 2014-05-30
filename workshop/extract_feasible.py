#!/usr/bin/python
def get_value(m, r):
	if m == '1':
		mod_value = 2;
	elif m == '2':
		mod_value = 4;
	else:
		mod_value = 0;
	return mod_value * float(r)

import os, os.path
fpath_prefix = './raw'
thres = 0.00001
feasible_list = list()
point_list = list()
d_point_list = list()
p_value = list()
p_mcs = list()
alpha_list = list()
for root, dirs, files in os.walk(fpath_prefix):
	for f in files:
		fullpath = os.path.join(root, f)
		#print fullpath
		f = open(fullpath, 'r')
		lines = f.readlines()
		d1 = lines[0].split()[1]
		d2 = lines[1].split()[1]
		d_point = (float(d1), float(d2))
		pl1 = lines[2].split()[1]
		pl2 = lines[3].split()[1]
		point = (float(pl1), float(pl2))
		mod1 = lines[4].split()[1]
		mod2 = lines[5].split()[1]
		rate1 = lines[6].split()[1]
		rate2 = lines[7].split()[1]
		mark = 0
		f_region = list()
		for l in lines:
			if mark == 1:
				#start collect
				words = l.split()
				if float(words[1]) < thres and float(words[2]) < thres:
					f_region.append(words[0])
			if l == 'alpha ber1 ber2\n':
				mark = 1
		feasible_list.append(f_region)
		if len(f_region) != 0:
			mid_fea = (float(f_region[0]) + float(f_region[-1])) / 2
			if point in point_list:
				#check if its better
				if get_value(mod1, rate1) + get_value(mod2, rate2) > p_value[point_list.index(point)]:
					p_value[point_list.index(point)] = get_value(mod1, rate1) + get_value(mod2, rate2)
					p_mcs[point_list.index(point)] = (mod1, mod2, rate1, rate2)
					alpha_list[point_list.index(point)] = mid_fea
			else:
				point_list.append(point)
				d_point_list.append(d_point)
				p_value.append(get_value(mod1, rate1) + get_value(mod2, rate2))
				p_mcs.append((mod1, mod2, rate1, rate2))
				alpha_list.append(mid_fea)
#output
output = open("ext_fea_out", "w")
print ''
print 'total: ', len(point_list)
for i in range(0, len(point_list)):
	print point_list[i], p_mcs[i], p_value[i]
	for p in point_list[i]:
		output.write(str(p) + " ")
	for mcs in p_mcs[i]:
		output.write(str(mcs) + " ")
	output.write(str(p_value[i]) + " ")
	output.write(str(alpha_list[i]) + " ")
	for d in d_point_list[i]:
		output.write(str(d) + " ")
	output.write("\n")
output.close()
