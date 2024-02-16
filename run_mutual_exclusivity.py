from network_generation import *
from network_statistics import *
import os

def main():
    xtcs = []
    for file in os.listdir('R1-15-closed'):
        if file.endswith('.xtc'):
            xtcs.append('R1-15-closed/'+file)
    xtcs.sort(key=lambda x: int(x.split('-')[1]))


    u = mda.Universe('R1-15-closed/R1-0-start-membrane-3JYC.pdb',*xtcs)

    mutual_exclusivity(u)

if __name__ == "__main__":
    main()
