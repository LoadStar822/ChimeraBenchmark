"""Re-evaluate saved quantification without rerunning read classification.

The count baseline is checked against the archived evaluation before native EM
profiles are scored against the identical full-rank truth and taxonomy.
"""
import argparse, csv, hashlib, json, shutil, sys
from pathlib import Path
from multiprocessing import Pool

CAP = {}

def worker(name):
    from chimera_bench.core import metrics as M
    from chimera_bench.tools.centrifuger_convert import convert_centrifuger_quant
    original_uni = M.compute_weighted_unifrac
    original_opal = M.compute_opal_profile_metrics
    def cap_uni(truth,pred,taxonomy):
        CAP.update(truth_full=truth,taxonomy=taxonomy)
        return original_uni(truth,pred,taxonomy)
    def cap_opal(truth,pred,ranks):
        CAP.update(truth_rank=truth,ranks=ranks)
        return original_opal(truth,pred,ranks)
    source=ARGS.runs/name
    target=ARGS.output/'runs'/name
    target.mkdir(parents=True,exist_ok=True)
    quant=ARGS.quant/name/'centrifuger_quant_reads.tsv'
    hashes={'quant':hashlib.sha256(quant.read_bytes()).hexdigest(), **{x:hashlib.sha256((source/x).read_bytes()).hexdigest() for x in ['meta.json','metrics.json']}}
    converter_sha = hashlib.sha256((ARGS.repo/'chimera_bench/tools/centrifuger_convert.py').read_bytes()).hexdigest()
    if (target/'DONE').exists():
        saved = json.loads((target/'metrics.json').read_text())['profile_provenance']
        if saved['sha256'] != hashes or saved['converter_sha256'] != converter_sha:
            raise ValueError(f'{name}: inputs or converter changed; use a new output directory')
        return name, 'cached'
    counts=target/'read_counts.cami.tsv'; native=target/'native_sequence.cami.tsv'
    convert_centrifuger_quant(quant_path=quant,out_path=counts,profile_method='read-counts')
    convert_centrifuger_quant(quant_path=quant,out_path=native)
    old=json.loads((source/'metrics.json').read_text()); ev=old['evaluation']; ds=dict(ev['dataset'])
    if 'truth_profile' in ds:ds['truth_profile']=ds['truth_profile'].replace('zymo_d6331.tsv','zymo_d6310.tsv')
    M.compute_weighted_unifrac=cap_uni
    M.compute_opal_profile_metrics=cap_opal
    try:
        baseline=M.evaluate_with_truth(ev['config'],ds,{'cami_profile_tsv':str(counts)},include_per_read=False)
    finally:
        M.compute_weighted_unifrac=original_uni
        M.compute_opal_profile_metrics=original_opal
    archived=json.loads((ARGS.quant/name/'eval.json').read_text())
    for k in ['l1_norm_species','weighted_unifrac','purity_species','completeness_species','l1_norm_genus','purity_genus','completeness_genus']:
        assert abs(baseline[k]-archived[k])<1e-9,(name,k,baseline[k],archived[k])
    pred_full=M.flatten_rank_profile(M.parse_cami_profile(native),CAP['taxonomy'])
    pred_rank,resolved=M.project_profile_to_ranks(pred_full,CAP['taxonomy'],CAP['ranks'])
    got=original_opal(CAP['truth_rank'],pred_rank,CAP['ranks'])
    got['weighted_unifrac']=original_uni(CAP['truth_full'],pred_full,CAP['taxonomy'])
    for rank,val in resolved.items():got['profile_resolved_fraction_'+rank]=val
    # Retain original per-read results; replace only profile evaluation fields.
    updated=dict(old); updated.update(baseline); updated.update(got)
    updated['evaluation']={**ev,'dataset':ds}
    updated['profile_method']='centrifuger-native-sequence-v1'
    updated['profile_provenance']={'original_run':str(source),'quant_path':str(quant),'sha256':hashes,'conversion':'species abundance * genomeSize, normalized over species; positive saved abundance; no additional threshold','converter_sha256':converter_sha}
    (target/'metrics.json').write_text(json.dumps(updated,indent=2)+'\n')
    shutil.copy2(source/'meta.json',target/'meta.json')
    (target/'read_counts_metrics.json').write_text(json.dumps(baseline,indent=2)+'\n')
    (target/'truth_full.json').write_text(json.dumps(CAP['truth_full'])+'\n')
    tv=M._normalize_profile({k:v for k,v in CAP['truth_rank']['species'].items() if isinstance(k,int)})
    pv=M._normalize_profile({k:v for k,v in pred_rank['species'].items() if isinstance(k,int)})
    with (target/'profile_species.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['taxid','truth','pred'])
        for t in sorted(tv.keys()|pv.keys()):w.writerow([t,tv.get(t,0),pv.get(t,0)])
    (target/'DONE').write_text('baseline validated; native profile evaluated\n')
    return name,got['l1_norm_species']

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--runs',type=Path,required=True);p.add_argument('--quant',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4);p.add_argument('--limit',type=int)
    ARGS=p.parse_args();sys.path.insert(0,str(ARGS.repo));names=sorted(d.name for d in ARGS.quant.iterdir() if (d/'centrifuger_quant_reads.tsv').is_file())
    if ARGS.limit:names=names[:ARGS.limit]
    with Pool(ARGS.workers) as pool:
        for i,result in enumerate(pool.imap_unordered(worker,names),1):print(i,*result,flush=True)
