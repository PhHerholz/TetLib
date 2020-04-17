import copy
import pandas as pd

def join_run_values(vdict, run_orbitpoints, rnames):
    val_dict = copy.deepcopy(vdict)
    for run_name in rnames:
        _, cellSize, rerat, facetSize, lloyd, perturb, exude, minorbitpoints, nflips, _ = run_name.split('_')
        n_entries = val_dict[run_name].shape[0]
        lloyd   = bool(int(lloyd))
        perturb = bool(int(perturb))
        exude   = bool(int(exude))
        nflips  = int(nflips)
        cellSize = float(cellSize)
        facetSize = float(facetSize)

        minorbitpoints = int(minorbitpoints)
        n_orbitpoints = len(run_orbitpoints[run_name]) - 1
        max_shift     = run_orbitpoints[run_name][' moved_by'].max()

        val_dict[run_name]['cellSize'] = [cellSize] * n_entries
        val_dict[run_name]['rerat']    = [rerat]    * n_entries
        val_dict[run_name]['lloyd']    = [lloyd] * n_entries
        val_dict[run_name]['perturb']  = [perturb] * n_entries
        val_dict[run_name]['exude']    = [exude] * n_entries
        val_dict[run_name]['nflips']   = [nflips] * n_entries
        val_dict[run_name]['minorbitpoints'] = [minorbitpoints] * n_entries
        val_dict[run_name]['n_orbitpoints'] = [n_orbitpoints] * n_entries
        val_dict[run_name]['max_shift'] = [max_shift] * n_entries

    return pd.concat([val_dict[run_name] for run_name in rnames], ignore_index=True)