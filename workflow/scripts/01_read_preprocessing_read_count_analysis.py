import pandas as pd
from plotnine import *

# read kneaddata read count data
kneaddata_read_counts = pd.read_csv(str(snakemake.input.kneaddata), sep='\t')

# if clumpify was run, incorporate this data as well
if snakemake.params.run_clumpify:
    clumpify_read_counts = pd.read_csv(str(snakemake.input.clumpify), sep='\t')
    read_counts = clumpify_read_counts.merge(kneaddata_read_counts, on='Sample', how='right')
    read_counts.rename(columns={'raw pair1': 'deduplicated pair1', 'raw pair2': 'deduplicated pair2', 'raw reads':'raw pair1'}, inplace=True)

    # melt read counts data into a table for plotting
    read_counts_melt = read_counts.melt(id_vars=['Sample'], value_vars=['raw pair1', 'deduplicated pair1', 'trimmed pair1', 'decontaminated hg37dec_v0.1 pair1', 'final pair1', 'final orphan1', 'final orphan2'])
    read_counts_melt['variable'] = pd.Categorical(read_counts_melt['variable'], categories=['raw pair1', 'deduplicated pair1', 'trimmed pair1', 'decontaminated hg37dec_v0.1 pair1', 'final pair1', 'final orphan1', 'final orphan2'], ordered = True)

else:
    # don't read clumpify data if set to false
    read_counts = kneaddata_read_counts

    # set order of categories when plotting
    read_counts_melt = read_counts.melt(id_vars=['Sample'], value_vars=['raw pair1','trimmed pair1', 'decontaminated hg37dec_v0.1 pair1', 'final pair1', 'final orphan1', 'final orphan2'])
    read_counts_melt['variable'] = pd.Categorical(read_counts_melt['variable'], categories=['raw pair1', 'trimmed pair1', 'decontaminated hg37dec_v0.1 pair1', 'final pair1', 'final orphan1', 'final orphan2'], ordered = True)

# save combined read counts file
read_counts.to_csv(str(snakemake.output.report), index=False, sep='\t')

# plot read counts data
read_counts_plot = (
    ggplot(read_counts_melt)
    + geom_boxplot(aes(x='variable', y='value'))
    + theme(figure_size=(16, 8)) 
    + labs(title='Read counts at each preprocessing step')
    + ylab("Read counts")
    + xlab("Preprocessing step")
)

# save figure to file
read_counts_plot.save(str(snakemake.output.figure), dpi=600)