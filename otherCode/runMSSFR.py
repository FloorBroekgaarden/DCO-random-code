import sys
import fileinput




import sys

def changeMetallicityRun(z, dirname):

    # replace all occurrences of 'sit' with 'SIT' and insert a line after the 5th
    for i, line in enumerate(fileinput.input('./masterFolder/pythonSubmit.py', inplace=1)):
        sys.stdout.write(line.replace('metallicity = 0.0142', 'metallicity = %s'%z))  



    # replace all occurrences of 'sit' with 'SIT' and insert a line after the 5th
    for i, line in enumerate(fileinput.input('compas_hpc_input.py', inplace=1)):
        sys.stdout.write(line.replace('/home/floorb/DATA/black_hole-neutron_star/', '/home/floorb/DATA/black_hole-neutron_star/' + str(dirname) +'/Z_%s'%str(z)))  # replace 'sit' and write





if __name__ == "__main__":
    z = float((sys.argv[1]))
    dirname = (sys.argv[2])
#    print('test')
    changeMetallicityRun(z,dirname)
