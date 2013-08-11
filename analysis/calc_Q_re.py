import os

n_mon = raw_input("Number of monomers? ")
density = raw_input("Solvent density? ")
a_cis = raw_input("Cis solvent interaction strength? ")
a_trans = raw_input("Trans solvent interaction strength? ")

data_path = '/Users/taylor/Documents/workspace/dpd_mc/data/v4.5/n' + n_mon + 'd' + density + 'c' + a_cis + 't' + a_trans

print data_path

try:
    os.chdir(data_path)
except:
    print "Could not open data folder"

Q_path = os.getcwd() + '/Q'
re_path = os.getcwd() + '/re'

re_z = [0] * 1001
count = [0] * 1001


for i in range(1,100):
    Q = []
    re = []

    Q_file = Q_path + '/Q' + str(i) + '.dat'
    re_file = re_path + '/re' + str(i) + '.dat'

    for line in open(Q_file):
        Q.append(float(line.split()[1]))

    for line in open(re_file):
        re.append(float(line.split()[4]))

    for j in range(1,len(Q)):
        index = int((Q[j] + 0.0005) / 0.001)
        count[index] += 1
        re_z[index] += abs(re[j])

    print "Finished " + re_file

file_out = open('Q_rez.dat', 'w')
for i in range(1,1000):
    if count[i] > 0:
        file_out.write("%f %f\n" % (i*0.001, re_z[i]/count[i]))

file_out.close()
    

