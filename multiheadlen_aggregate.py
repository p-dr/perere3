import utils as u
import aggregate_data

## WARNING: aggregate_data.outpath is not reset in the end
## WARNING: headlengths are wrong in the final aggregated at outagg
for count_dir in u.pardir.glob('*_counted_reads'):
    prefix = count_dir.stem.split('_')[0]
    u.log(f'Aggregating {prefix}')

    outagg = aggregate_data.outpath 
    outagg = outagg.parent/(prefix + '_' + outagg.name.split('bp_')[-1])
    aggregate_data.outpath = outagg

    aggregate_data.main()
