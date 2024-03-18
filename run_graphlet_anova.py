from graphlets import *
import warnings
from scipy.stats import f_oneway
import statsmodels.api as sm
from statsmodels.formula.api import ols





def main():
    ldirs = [f'../orca/output/{r}-{c}-closed' for r in ['R1','R2','R3'] for c in ['15','30']]
    # load data from uniform slices into df from multiple replicates
    gdds = []
    for d in ldirs:
        for f in os.listdir(d):
            if f == '.DS_Store': continue
            if '-0-' in f: continue # skip the one with the burn in time
            r = f.split('-')[0]
            gdd = graphlet_degree_distribution(os.path.join(d, f))

            if '30' in f: chol = 30
            else: chol = 15
            if 'closed' in f: state = 'closed'
            else: state = 'open'
            gdds.append([f,chol,state,r] + list(gdd))

    gcols = ['g' + str(graphlet) for graphlet in range(73)]
    df = pd.DataFrame(gdds,columns=["name", "chol", "state",'replicate'] +gcols)
    # ANOVA test for each coordinate of vector.
    # don't want to use MANOVA because want to evaluate each graphlet separately
    
    keys = []
    tables = []
    for g in range(73):
        model = ols(f'g{g} ~ chol', data=df).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)

        keys.append(g)
        tables.append(anova_table)

    df_anova = pd.concat(tables, keys=keys, axis=0)
    df_anova.to_csv('df_anova.csv',index=False)

    # do for whole network and for nodes.



if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    main()