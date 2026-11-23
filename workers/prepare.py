"""Verify the two source FASTQs and prepare all records in one matched library."""
import csv
import json
import re
import sys
from pathlib import Path

root = Path(__file__).resolve().parent
agent = root.parent
root = Path.cwd() if (Path.cwd()/"pair.json").is_file() else root
pair = json.loads((root/"pair.json").read_text()) if (root/"pair.json").is_file() else {"donor":"D3734","replicate":"T1"}
sys.path.insert(0, str(agent))
from absolve.settings import MANIFEST, REFERENCES
from tools.pilot import FROZEN, PUBLIC, one, public_library, table
sys.path.insert(0, str(FROZEN/'src'))
import prepare_crispresso as prep
import barcode_core as bc

condition = pair.get('conditions', ('NoEP','HiFi'))[int(sys.argv[1])]
donor, replicate = pair['donor'], pair['replicate']
with MANIFEST.open() as handle:
    row = one(r for r in csv.DictReader(handle,delimiter='\t') if r['sample_name']==f'{condition}_pool_{donor}_{replicate}' and r['status']=='matched')
cas_type = {'NoEP':'noEP','HiFi':'HiFiCas9','WT':'WTCas9'}[condition]
meta = one(r for r in table(PUBLIC/'NovaSeq3_sample_info.csv')
           if (r['donor_id'],r['replicate'],r['cas_type'])==(donor,replicate,cas_type))
sample = row['sample_name'] if condition == 'WT' else meta['LV_ID']
archive_library = public_library(row)
if archive_library and archive_library != sample:
    raise ValueError('Public run and reference sample table disagree: ' + row['sample_name'])
work=root/condition
work.mkdir()
record={'pool_id':sample}
for mate in ('r1','r2'):
    path=row[mate]
    name = Path(path).name
    if condition == 'WT':
        assert re.fullmatch(r'Nola_EP031022_\d+_LVpool_wtCas9_'+donor+'_'+replicate+'_R'+mate[-1]+r'\.fastq\.gz', name)
    elif archive_library is None:
        assert name.startswith(sample+'_CKDL')
    assert bc.sha256(path)==row['sha256_'+mate]
    record[mate]=path
    record[mate+'_identity']=bc.identity(path)
manifest=work/'inputs.json'
manifest.write_text(json.dumps([record],indent=2))
config=FROZEN/'config'
policy=json.loads((config/'barcode_policy.json').read_text())
policy['structures']={'original':''}
policy_file=work/'original_policy.json'
policy_file.write_text(json.dumps(policy))
if condition == 'WT':
    from wt_layout import WTRecognizer
    # Only this WT worker process uses the repository's fixed-region rule.
    bc.Recognizer = WTRecognizer
stats=prep.prepare(manifest,sample,policy_file,config/'crispresso_prep_policy.json',REFERENCES/'targets.json',work/'prepared',None)
assert stats['scope']=='full' and stats['counts_conserved']
print(condition,stats['pairs'],stats['prepared_pairs'],stats['targets_with_reads'])
