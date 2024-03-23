from dynamic_graphlets import *

def main():
    winlen = 50

    for r in ['R1','R2','R3']:
        print(f'working on {r}')
        # make dynamic graphlet input
        for i in ['30','15']:
            print(f'making graph inputs for {r}-{i}')
            xtcs = []
            for file in os.listdir(f'{r}-{i}-closed'):
                if file.endswith('.xtc'):
                    xtcs.append(f'{r}-{i}-closed/'+file)
            xtcs.sort(key=lambda x: int(x.split('-')[1]))
            u = mda.Universe(f'{r}-{i}-closed/{r}-0-start-membrane-3JYC.pdb',*xtcs)
            starts = np.linspace(0,len(u.trajectory)-winlen,100).astype(np.int32)
            basename = f'../dynamic_graphlets/input/{r}-{i}-closed'
            for j in tqdm.tqdm(starts):
                fn = f'{basename}-{j}-{j+winlen}.in'
                output_temporal_graph(u,fn,s=j,d=j+winlen)



if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()