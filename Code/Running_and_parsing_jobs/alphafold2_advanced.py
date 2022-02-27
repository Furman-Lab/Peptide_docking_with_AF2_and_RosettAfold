# -*- coding: utf-8 -*-
"""AlphaFold2 running script based on AlphaFold2_advanced.ipynb from ColabFold repository

Original file is located at
    https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/beta/AlphaFold2_advanced.ipynb
See [ColabFold](https://github.com/sokrypton/ColabFold/) for other related notebooks by Sergey Ovchinnikov

"""
import argparse
import os
os.environ['TF_FORCE_UNIFIED_MEMORY']='1'
os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION']='2.0'
os.environ['TF_FORCE_GPU_ALLOW_GROWTH ']='true'

#### Read arguments ####
parser = argparse.ArgumentParser(description='Call AlphaFold with advanced parameters')
group_input = parser.add_argument_group('Input')
group_output = parser.add_argument_group('Output')
group_msa = parser.add_argument_group('Control MSA')
group_model = parser.add_argument_group('Control model parameters')
group_relax = parser.add_argument_group('Relax parameters')

group_input.add_argument('-s', dest='sequence', action='store', help='The sequence(s) used for prediction. Use / to specify intra-protein chainbreaks (for trimming regions within protein). Use : to specify inter-protein chainbreaks (for modeling protein-protein hetero-complexes).', default=None)
group_input.add_argument('--fasta', dest='fasta', action='store', help='The sequences for prediction, one per', default=None)
group_input.add_argument('--num_models', dest='num_models', type=int, default=5, help='Number of models.')
group_input.add_argument('--homooligomer', dest='homooligomer', type=str, default='1', help='Number of times to repeat the CONCATENATED sequence')

group_output.add_argument('--prefix', dest='prefix', help='Prefix for every file, e.g. complex name')
group_output.add_argument('--working_dir', dest='working_dir', default='.', help='Working directory')
group_output.add_argument('--dpi', dest='dpi', type=int, default=100, help='DPI of output figures')
group_output.add_argument('--save_pae_json', dest='save_pae_json', action='store_true', default=False, help='')
group_output.add_argument('--delete_files', dest='delete_files', action='store_true', default=False, help='Delete old files before prediction')

group_model.add_argument('--max_recycles', dest='max_recycles', type=int, default=3, help='Controls the maximum number of times the structure is fed back into the neural network for refinement.')
group_model.add_argument('--use_ptm', dest='use_ptm', action='store_true', default=False,
          help="Uses Deepmind's `ptm` finetuned model parameters to get PAE per structure. Disable to use the original model params.")
group_model.add_argument('--is_training', dest='is_training', action='store_true', default=False,
          help="enables the stochastic part of the model (dropout), when coupled with `num_samples` can be used to sample a diverse set of structures")
group_model.add_argument('--tol', dest='tol', type=int, default=0, help='Tolerance for deciding when to stop (CA-RMS between recycles)')
group_model.add_argument('--num_ensemble', dest='num_ensemble', type=int, default=1,
          help='The trunk of the network is run multiple times with different random choices for the MSA cluster centres. (`1`=`default`, `8`=`casp14 setting`)')
group_model.add_argument('--num_samples', dest='num_samples', type=int, default=1,
          help='Sets number of random_seeds to iterate through for each model.')
group_model.add_argument('--use_turbo', dest='use_turbo', action='store_true', default=True,
          help='Introduces a few modifications (compile once, swap params, adjust max_msa) to speedup and reduce memory requirements. Disable for default behavior.')

group_msa.add_argument('--subsample_msa', dest='subsample_msa', action='store_true', help='subsample msa to avoid OOM error')
group_msa.add_argument('--just_msa', dest='just_msa', action='store_true', help='create MSA and exit')
group_msa.add_argument('--no_use_env', dest='no_use_env', action='store_true', help='use environmental sequences?')
group_msa.add_argument('--pair_msa', dest='pair_msa', action='store_true', help='pair msa for prokaryotic sequences')
group_msa.add_argument('--msa_method', dest='msa_method', type=str, default="mmseqs2", help='MSA method. [mmseqs2, jackhammer, single_sequence, custom_a3m, precomputed]')
group_msa.add_argument('--custom_a3m', dest='custom_a3m', type=str, help='In case msa_method=custom_a3m, this option is required')
group_msa.add_argument('--rank_by', dest='rank_by', type=str, default='pLDDT', help='specify metric to use for ranking models (For protein-protein complexes, we recommend pTMscore). [pLDDT, pTMscore]')
group_msa.add_argument('--max_msa', dest='max_msa', type=str, default='512:1024', help='defines: `max_msa_clusters:max_extra_msa` number of sequences to use. (Lowering will reduce GPU requirements, but may result in poor model quality.')
group_msa.add_argument('--precomputed', dest='precomputed', type=str,
          help='In case msa_method=precomputed, this option is required. If you have previously run this notebook and saved the results, you can skip MSA generating step by providing the previously generated  `prediction/msa.npz')
group_msa.add_argument('--cov', dest='cov', type=int, default=0, help="Filter to remove sequences that don't cover at least `cov` %% of query. (Set to `0` to disable all filtering.)")

group_relax.add_argument('--use_amber_relax', dest='use_amber_relax', action='store_true', default=False, help='Use amber relaxation - NOT YET SUPPORTED')
group_relax.add_argument('--relax_all', dest='relax_all', action='store_true', default=False,
          help='Amber-relax all models. Disable to only relax the top ranked model. (Note: no models are relaxed if `use_amber_relax` is disabled. - NOT YET SUPPORTED')
args = parser.parse_args()
print(args)

if(args.msa_method=='custom_a3m' and args.custom_a3m is None):
  print('If custom_a3m is set, a file for option --custom_a3m must be provided. See -h for more information. Exiting.')
  exit()
if(args.fasta==None and args.sequence is None):
  print('Provide sequences for prediction either as sequences or in FASTA format. Exiting.')
  exit()
if(not args.fasta==None and not args.sequence is None):
  print('Provide sequences for prediction EITHER as sequences or in FASTA format, but not both. Exiting.')
  exit()
if(args.msa_method=='precomputed' and args.precomputed is None):
  print('If precomputed is set, a file for option --precomputed must be provided. See -h for more information. Exiting.')
  exit()
if(args.msa_method=='custom_a3m' and args.custom_a3m is None):
  print('If custom_a3m is set, a file for option --custom_a3m must be provided. See -h for more information. Exiting.')
  exit()

rank_by = args.rank_by
max_msa = args.max_msa
max_msa_clusters, max_extra_msa = [int(x) for x in max_msa.split(":")]
just_msa = args.just_msa

use_amber_relax = args.use_amber_relax # NOT SUPPORTED, need more installation
use_ptm = args.use_ptm
max_recycles = args.max_recycles
num_models = args.num_models
tol = args.tol
num_ensemble = args.num_ensemble
num_samples = args.num_samples
pair_msa = args.pair_msa
subsample_msa = args.subsample_msa

use_turbo = args.use_turbo
relax_all = args.relax_all
save_pae_json = args.save_pae_json
dpi = args.dpi
cov = args.cov
msa_method = args.msa_method
homooligomer = args.homooligomer
is_training = args.is_training

if args.no_use_env == True:
	use_env = False
else:
	use_env = True

MIN_SEQUENCE_LENGTH = 16
MAX_SEQUENCE_LENGTH = 2500

##############################################################x
# prepare sequence
import re
sequence = args.sequence
sequence = re.sub("[^A-Z:/]", "", sequence.upper())
sequence = re.sub(":+",":",sequence)
sequence = re.sub("/+","/",sequence)

# define number of copies
homooligomer =  "1" #@param {type:"string"}
if len(homooligomer) == 0: homooligomer = "1"
homooligomer = re.sub("[^0-9:]", "", homooligomer)
homooligomers = [int(h) for h in homooligomer.split(":")]


ori_sequence = sequence
sequence = sequence.replace("/","").replace(":","")
seqs = ori_sequence.replace("/","").split(":")

# prepare homo-oligomeric sequence
if len(seqs) != len(homooligomers):
  if len(homooligomers) == 1:
    homooligomers = [homooligomers[0]] * len(seqs)
    homooligomer = ":".join([str(h) for h in homooligomers])
  else:
    while len(seqs) > len(homooligomers):
      homooligomers.append(1)
    homooligomers = homooligomers[:len(seqs)]
    homooligomer = ":".join([str(h) for h in homooligomers])
    print("WARNING: Mismatch between number of breaks ':' in 'sequence' and 'homooligomer' definition")

full_sequence = "".join([s*h for s,h in zip(seqs,homooligomers)])

# prediction directory
output_dir = args.working_dir # 'prediction_' + cf.get_hash(full_sequence)[:5]
os.makedirs(output_dir, exist_ok=True)
print(f"working directory: {output_dir}")

# print out params
print(f"homooligomer: '{homooligomer}'")
print(f"total_length: '{len(full_sequence)}'")
print(f"working_directory: '{output_dir}'")

########################
# --- Python imports ---
import tensorflow as tf
import jax
from IPython.utils import io
import subprocess
print(os.getcwd())
import colabfold as cf
import sys
import pickle
if use_amber_relax:
  sys.path.append('/opt/conda/lib/python3.7/site-packages')

from urllib import request
from concurrent import futures
import json
from matplotlib import gridspec
import matplotlib.pyplot as plt
import numpy as np

from alphafold.model import model
from alphafold.model import config
from alphafold.model import data

from alphafold.data import parsers
from alphafold.data import pipeline
from alphafold.data.tools import jackhmmer

from alphafold.common import protein

#############
gpus = tf.config.list_physical_devices('GPU')
if gpus:
  try:
    # Currently, memory growth needs to be the same across GPUs
    for gpu in gpus:
      print(gpu)
      tf.config.experimental.set_memory_growth(gpu, True)
    logical_gpus = tf.config.list_logical_devices('GPU')
    print(len(gpus), "Physical GPUs,", len(logical_gpus), "Logical GPUs")
  except RuntimeError as e:
    # Memory growth must be set before GPUs have been initialized
    print(e)
############

if use_amber_relax:
  from alphafold.relax import relax
  from alphafold.relax import utils

def run_jackhmmer(sequence, prefix):
  pickled_msa_path = f"{prefix}.jackhmmer.pickle"
  if os.path.isfile(pickled_msa_path):
    msas_dict = pickle.load(open(pickled_msa_path,"rb"))
    msas, deletion_matrices = (msas_dict[k] for k in ['msas', 'deletion_matrices'])
    full_msa = []
    for msa in msas:
      full_msa += msa
  else:
    # --- Find the closest source ---
    test_url_pattern = 'https://storage.googleapis.com/alphafold-colab{:s}/latest/uniref90_2021_03.fasta.1'
    ex = futures.ThreadPoolExecutor(3)
    def fetch(source):
      request.urlretrieve(test_url_pattern.format(source))
      return source
    fs = [ex.submit(fetch, source) for source in ['', '-europe', '-asia']]
    source = None
    for f in futures.as_completed(fs):
      source = f.result()
      ex.shutdown()
      break

    jackhmmer_binary_path = '/vol/ek/share/bin/hmmer-3.1b1/bin/jackhmmer'
    dbs = []

    num_jackhmmer_chunks = {'uniref90': 59, 'smallbfd': 17, 'mgnify': 71}
    total_jackhmmer_chunks = sum(num_jackhmmer_chunks.values())

  print('Searching uniref90')
  jackhmmer_uniref90_runner = jackhmmer.Jackhmmer(
    binary_path=jackhmmer_binary_path,
    database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/uniref90_2021_03.fasta',
    get_tblout=True,
    num_streamed_chunks=num_jackhmmer_chunks['uniref90'],
    streaming_callback=jackhmmer_chunk_callback,
    z_value=135301051)
  dbs.append(('uniref90', jackhmmer_uniref90_runner.query('target.fasta')))

  print('Searching smallbfd')
  jackhmmer_smallbfd_runner = jackhmmer.Jackhmmer(
    binary_path=jackhmmer_binary_path,
    database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/bfd-first_non_consensus_sequences.fasta',
    get_tblout=True,
    num_streamed_chunks=num_jackhmmer_chunks['smallbfd'],
    streaming_callback=jackhmmer_chunk_callback,
    z_value=65984053)
  dbs.append(('smallbfd', jackhmmer_smallbfd_runner.query('target.fasta')))

  print('Searching mgnify')
  jackhmmer_mgnify_runner = jackhmmer.Jackhmmer(
    binary_path=jackhmmer_binary_path,
    database_path=f'https://storage.googleapis.com/alphafold-colab{source}/latest/mgy_clusters_2019_05.fasta',
    get_tblout=True,
    num_streamed_chunks=num_jackhmmer_chunks['mgnify'],
    streaming_callback=jackhmmer_chunk_callback,
    z_value=304820129)
  dbs.append(('mgnify', jackhmmer_mgnify_runner.query('target.fasta')))

  # --- Extract the MSAs and visualize ---
  # Extract the MSAs from the Stockholm files.
  # NB: deduplication happens later in pipeline.make_msa_features.

  mgnify_max_hits = 501
  msas = []
  deletion_matrices = []
  for db_name, db_results in dbs:
    unsorted_results = []
    for i, result in enumerate(db_results):
      msa, deletion_matrix, target_names = parsers.parse_stockholm(result['sto'])
      e_values_dict = parsers.parse_e_values_from_tblout(result['tbl'])
      e_values = [e_values_dict[t.split('/')[0]] for t in target_names]
      zipped_results = zip(msa, deletion_matrix, target_names, e_values)
      if i != 0:
        # Only take query from the first chunk
        zipped_results = [x for x in zipped_results if x[2] != 'query']
      unsorted_results.extend(zipped_results)
    sorted_by_evalue = sorted(unsorted_results, key=lambda x: x[3])
    db_msas, db_deletion_matrices, _, _ = zip(*sorted_by_evalue)
    if db_msas:
      if db_name == 'mgnify':
        db_msas = db_msas[:mgnify_max_hits]
        db_deletion_matrices = db_deletion_matrices[:mgnify_max_hits]
      msas.append(db_msas)
      deletion_matrices.append(db_deletion_matrices)
      msa_size = len(set(db_msas))
      print(f'{msa_size} Sequences Found in {db_name}')

    pickle.dump({"msas":msas,"deletion_matrices":deletion_matrices},
                open(pickled_msa_path,"wb"))
  return msas, deletion_matrices


aatypes = set('ACDEFGHIKLMNPQRSTVWY')  # 20 standard aatypes
if not set(full_sequence).issubset(aatypes):
  raise Exception(f'Input sequence contains non-amino acid letters: {set(sequence) - aatypes}. AlphaFold only supports 20 standard amino acids as inputs.')
if len(full_sequence) < MIN_SEQUENCE_LENGTH:
  raise Exception(f'Input sequence is too short: {len(full_sequence)} amino acids, while the minimum is {MIN_SEQUENCE_LENGTH}')
if len(full_sequence) > MAX_SEQUENCE_LENGTH:
  raise Exception(f'Input sequence is too long: {len(full_sequence)} amino acids, while the maximum is {MAX_SEQUENCE_LENGTH}. Please use the full AlphaFold system for long sequences.')

if len(full_sequence) > 1400:
  print(f"WARNING: For a typical Google-Colab-GPU (16G) session, the max total length is ~1300 residues. You are at {len(full_sequence)}! Run Alphafold may crash.")

# tmp directory
prefix = cf.get_hash(sequence)
os.makedirs('tmp', exist_ok=True)
prefix = os.path.join('tmp',prefix)

# --- Search against genetic databases ---
with open('target.fasta', 'wt') as f:
  f.write(f'>query\n{sequence}')

# Run the search against chunks of genetic databases (since the genetic
# databases don't fit in Colab ramdisk).

if msa_method == "precomputed":
  print("use precomputed pickled msa from previous run")
  with open(args.precomputed, "rb") as msa_pickle:
    msas_dict = pickle.load(msa_pickle)
  msas, deletion_matrices = (msas_dict[k] for k in ['msas', 'deletion_matrices'])

elif msa_method == "single_sequence":
  msas = [[sequence]]
  deletion_matrices = [[[0] * len(sequence)]]

elif msa_method == "custom_a3m":
  print("use custom a3m")
  a3m_file = open(args.custom_a3m, "r")
  a3m_content = f.read()
  lines = a3m_content.splitlines()
  a3m_lines = []
  for line in lines:
    line = line.replace("\x00", "")
    if len(line) > 0 and not line.startswith('#'):
      a3m_lines.append(line)
  msa, deletion_matrix = parsers.parse_a3m("\n".join(a3m_lines))
  msas, deletion_matrices = [msa], [deletion_matrix]

  if len(msas[0][0]) != len(sequence):
    print("ERROR: the length of msa does not match input sequence")

else:
  os.makedirs('tmp', exist_ok=True)
  seqs = ori_sequence.replace('/', '').split(':')

  _blank_seq = ["-" * len(seq) for seq in seqs]
  _blank_mtx = [[0] * len(seq) for seq in seqs]


  def _pad(ns, vals, mode):
    if mode == "seq": _blank = _blank_seq.copy()
    if mode == "mtx": _blank = _blank_mtx.copy()
    if isinstance(ns, list):
      for n, val in zip(ns, vals): _blank[n] = val
    else:
      _blank[ns] = vals
    if mode == "seq": return "".join(_blank)
    if mode == "mtx": return sum(_blank, [])


  # gather msas
  msas, deletion_matrices = [], []
  if msa_method == "mmseqs2":
    prefix = cf.get_hash("".join(seqs))
    prefix = os.path.join('tmp', prefix)
    print(f"running mmseqs2")
    A3M_LINES = cf.run_mmseqs2(seqs, prefix, filter=True)

  for n, seq in enumerate(seqs):
    # tmp directory
    prefix = cf.get_hash(seq)
    prefix = os.path.join('tmp', prefix)

    if msa_method == "mmseqs2":
      # run mmseqs2
      a3m_lines = A3M_LINES[n]
      msa, mtx = parsers.parse_a3m(a3m_lines)
      msas_, mtxs_ = [msa], [mtx]

    elif msa_method == "jackhmmer":
      print(f"running jackhmmer on seq_{n}")
      # run jackhmmer
      msas_, mtxs_, names_ = run_jackhmmer(seq, prefix)

    # pad sequences
    for msa_, mtx_ in zip(msas_, mtxs_):
      msa, mtx = [sequence], [[0] * len(sequence)]
      for s, m in zip(msa_, mtx_):
        msa.append(_pad(n, s, "seq"))
        mtx.append(_pad(n, m, "mtx"))

      msas.append(msa)
      deletion_matrices.append(mtx)

  if pair_msa and len(seqs) > 1:
    print("attempting to pair some sequences...")

    if msa_method == "mmseqs2":
      prefix = cf.get_hash("".join(seq))
      prefix = os.path.join('tmp', prefix)
      print(f"running mmseqs2_noenv_nofilter on all seqs")
      A3M_LINES = cf.run_mmseqs2(seqs, prefix, use_env=False, filter=False)

    _data = []
    for a in range(len(seqs)):
      _seq = seqs[a]
      _prefix = os.path.join('tmp', cf.get_hash(_seq))

      if msa_method == "mmseqs2":
        a3m_lines = A3M_LINES[a]
        _msa, _mtx, _lab = pairmsa.parse_a3m(a3m_lines)

      elif msa_method == "jackhmmer":
        _msas, _mtxs, _names = run_jackhmmer(_seq, _prefix)
        _msa, _mtx, _lab = pairmsa.get_uni_jackhmmer(_msas[0], _mtxs[0], _names[0])

      if len(_msa) > 1:
        _data.append(pairmsa.hash_it(_msa, _lab, _mtx, call_uniprot=False))
      else:
        _data.append(None)

    for a in range(len(seqs)):
      if _data[a] is not None:
        for b in range(a + 1, len(seqs)):
          if _data[b] is not None:
            _seq_a, _seq_b, _mtx_a, _mtx_b = pairmsa.stitch(_data[a], _data[b])
            print(f"attempting to pair seq_{a} and seq_{b}... FOUND: {len(_seq_a)}")
            if len(_seq_a) > 0:
              msa, mtx = [sequence], [[0] * len(sequence)]
              for s_a, s_b, m_a, m_b in zip(_seq_a, _seq_b, _mtx_a, _mtx_b):
                msa.append(_pad([a, b], [s_a, s_b], "seq"))
                mtx.append(_pad([a, b], [m_a, m_b], "mtx"))
              msas.append(msa)
              deletion_matrices.append(mtx)

# save MSA as pickle
with open(os.path.join(output_dir,"msa.pickle"), "wb") as output_file:
    pickle.dump({"msas":msas,"deletion_matrices":deletion_matrices}, output_file)


if just_msa:
    print('MSA created, exiting...')
    exit()

if msa_method != "single_sequence" and cov > 0:
  # filter sequences that don't cover at least %
  msas, deletion_matrices = cf.cov_filter(msas, deletion_matrices, cov)

full_msa = []
for msa in msas: full_msa += msa

# deduplicate
deduped_full_msa = list(dict.fromkeys(full_msa))
total_msa_size = len(deduped_full_msa)
if msa_method == "mmseqs2":
  print(f'\n{total_msa_size} Sequences Found in Total (after filtering)\n')
else:
  print(f'\n{total_msa_size} Sequences Found in Total\n')

msa_arr = np.array([list(seq) for seq in deduped_full_msa])
num_alignments, num_res = msa_arr.shape

if num_alignments > 1:
  plt.figure(figsize=(8,5),dpi=100)
  plt.title("Sequence coverage")
  seqid = (np.array(list(sequence)) == msa_arr).mean(-1)
  seqid_sort = seqid.argsort() #[::-1]
  non_gaps = (msa_arr != "-").astype(float)
  non_gaps[non_gaps == 0] = np.nan
  plt.imshow(non_gaps[seqid_sort]*seqid[seqid_sort,None],
            interpolation='nearest', aspect='auto',
            cmap="rainbow_r", vmin=0, vmax=1, origin='lower',
            extent=(0, msa_arr.shape[1], 0, msa_arr.shape[0]))
  plt.plot((msa_arr != "-").sum(0), color='black')
  plt.xlim(0,msa_arr.shape[1])
  plt.ylim(0,msa_arr.shape[0])
  plt.colorbar(label="Sequence identity to query",)
  plt.xlabel("Positions")
  plt.ylabel("Sequences")
  plt.savefig(os.path.join(output_dir,"msa_coverage.png"), bbox_inches = 'tight', dpi=200)
  plt.show()

from string import ascii_uppercase

if use_ptm == False and rank_by == "pTMscore":
  print("WARNING: models will be ranked by pLDDT, 'use_ptm' is needed to compute pTMscore")
  rank_by = "pLDDT"

#############################
# delete old files
#############################
if args.delete_files:
	for f in os.listdir(output_dir):
	  if "rank_" in f:
	    os.remove(os.path.join(output_dir, f))

#############################
# homooligomerize
#############################
lengths = [len(seq) for seq in seqs]
msas_mod, deletion_matrices_mod = cf.homooligomerize_heterooligomer(msas, deletion_matrices,
                                                                    lengths, homooligomers)
#############################
# define input features
#############################
def _placeholder_template_feats(num_templates_, num_res_):
  return {
      'template_aatype': np.zeros([num_templates_, num_res_, 22], np.float32),
      'template_all_atom_masks': np.zeros([num_templates_, num_res_, 37, 3], np.float32),
      'template_all_atom_positions': np.zeros([num_templates_, num_res_, 37], np.float32),
      'template_domain_names': np.zeros([num_templates_], np.float32),
      'template_sum_probs': np.zeros([num_templates_], np.float32),
  }

num_res = len(full_sequence)
feature_dict = {}
feature_dict.update(pipeline.make_sequence_features(full_sequence, 'test', num_res))
feature_dict.update(pipeline.make_msa_features(msas_mod, deletion_matrices=deletion_matrices_mod))
if not use_turbo:
  feature_dict.update(_placeholder_template_feats(0, num_res))

def subsample_msa(F, N=10000, random_seed=0):
  '''subsample msa to avoid running out of memory'''
  M = len(F["msa"])
  if N is not None and M > N:
    np.random.seed(random_seed)
    idx = np.append(0,np.random.permutation(np.arange(1,M)))[:N]
    F_ = {}
    F_["msa"] = F["msa"][idx]
    F_["deletion_matrix_int"] = F["deletion_matrix_int"][idx]
    F_["num_alignments"] = np.full_like(F["num_alignments"],N)
    for k in ['aatype', 'between_segment_residues',
              'domain_name', 'residue_index',
              'seq_length', 'sequence']:
              F_[k] = F[k]
    return F_
  else:
    return F

################################
# set chain breaks
################################
Ls = []
for seq, h in zip(ori_sequence.split(":"),homooligomers):
  Ls += [len(s) for s in seq.split("/")] * h

Ls_plot = sum([[len(seq)]*h for seq, h in zip(seqs, homooligomers)], [])
feature_dict['residue_index'] = cf.chain_break(feature_dict['residue_index'], Ls)

###########################
# run alphafold
###########################
def parse_results(prediction_result, processed_feature_dict):
  b_factors = prediction_result['plddt'][:,None] * prediction_result['structure_module']['final_atom_mask']
  out = {"unrelaxed_protein": protein.from_prediction(processed_feature_dict, prediction_result, b_factors=b_factors),
         "plddt": prediction_result['plddt'],
         "pLDDT": prediction_result['plddt'].mean(),
         "dists": prediction_result["distogram"]["bin_edges"][prediction_result["distogram"]["logits"].argmax(-1)],
         "adj": jax.nn.softmax(prediction_result["distogram"]["logits"])[:,:,prediction_result["distogram"]["bin_edges"] < 8].sum(-1)}
  if "ptm" in prediction_result:
    out.update({"pae": prediction_result['predicted_aligned_error'],
                "pTMscore": prediction_result['ptm']})
  return out

model_names = ['model_1', 'model_2', 'model_3', 'model_4', 'model_5'][:num_models]
total = len(model_names) * num_samples
if use_amber_relax:
  if relax_all: total += total
  else: total += 1

#######################################################################
# precompile model and recompile only if length changes
if use_turbo:
  name = "model_5_ptm" if use_ptm else "model_5"
  N = len(feature_dict["msa"])
  L = len(feature_dict["residue_index"])
  compiled = (N, L, use_ptm, max_recycles, tol, num_ensemble, max_msa, is_training)
  if "COMPILED" in dir():
    if COMPILED != compiled: recompile = True
  else:
    recompile = True
  if recompile:
    cf.clear_mem("gpu")
    cfg = config.model_config(name)
    cfg.data.common.max_extra_msa = min(N, max_extra_msa)
    cfg.data.eval.max_msa_clusters = min(N, max_msa_clusters)
    cfg.data.common.num_recycle = max_recycles
    cfg.model.num_recycle = max_recycles
    cfg.model.recycle_tol = tol
    cfg.data.eval.num_ensemble = num_ensemble

    params = data.get_model_haiku_params(name, '/vol/ek/share/peptide_docking_with_afold2_and_rosettAfold/scripts/alphafold/data')
    model_runner = model.RunModel(cfg, params, is_training=is_training)
    COMPILED = compiled
    recompile = False
else:
  cf.clear_mem("gpu")
  recompile = True

# cleanup
if "outs" in dir(): del outs
outs = {}
cf.clear_mem("cpu")

#######################################################################
for num, model_name in enumerate(model_names):  # for each model
  name = model_name + "_ptm" if use_ptm else model_name

  # setup model and/or params
  params = data.get_model_haiku_params(name, '/vol/ek/share/peptide_docking_with_afold2_and_rosettAfold/scripts/alphafold/data')
  if use_turbo:
    for k in model_runner.params.keys():
      model_runner.params[k] = params[k]
  else:
    cfg = config.model_config(name)
    cfg.data.common.num_recycle = cfg.model.num_recycle = max_recycles
    cfg.model.recycle_tol = tol
    cfg.data.eval.num_ensemble = num_ensemble
    model_runner = model.RunModel(cfg, params, is_training=is_training)

  for seed in range(num_samples):  # for each seed -
    # predict
    key = f"{name}_seed_{seed}"

    if subsample_msa:
      print('Subsampling MSA')
      sampled_feats_dict = subsample_msa(feature_dict, random_seed=seed)
      processed_feature_dict = model_runner.process_features(sampled_feats_dict, random_seed=seed)
    else:
      processed_feature_dict = model_runner.process_features(feature_dict, random_seed=seed)

    # processed_feature_dict = model_runner.process_features(feature_dict, random_seed=seed)
    prediction_result, (r, t) = cf.to(model_runner.predict(processed_feature_dict, random_seed=seed), "cpu")
    outs[key] = parse_results(prediction_result, processed_feature_dict)

    # report
    line = f"{key} recycles:{r} tol:{t:.2f} pLDDT:{outs[key]['pLDDT']:.2f}"
    if use_ptm: line += f" pTMscore:{outs[key]['pTMscore']:.2f}"
    print(line)

    # cleanup
    del processed_feature_dict, prediction_result

  if use_turbo:
    del params
  else:
    del params, model_runner, cfg
    cf.clear_mem("gpu")

# delete old files
    for f in os.listdir(output_dir):
      if "rank_" in f:
        os.remove(os.path.join(output_dir, f))

# Find the best model according to the mean rank_by
model_rank = list(outs.keys())
model_rank = [model_rank[i] for i in np.argsort([outs[x][rank_by] for x in model_rank])[::-1]]

# Write out the prediction
for n,key in enumerate(model_rank):
  prefix = f"{args.prefix}_rank_{n+1}_{key}_recycle_{max_recycles}"
  pred_output_path = os.path.join(output_dir,f'{prefix}_unrelaxed.pdb')

  pdb_lines = protein.to_pdb(outs[key]["unrelaxed_protein"])
  with open(pred_output_path, 'w') as f:
    f.write(pdb_lines)
  if use_amber_relax:
    print(f'AMBER relaxation')
    if relax_all or n == 0:
      amber_relaxer = relax.AmberRelaxation(
          max_iterations=0,
          tolerance=2.39,
          stiffness=10.0,
          exclude_residues=[],
          max_outer_iterations=20)
      relaxed_pdb_lines, _, _ = amber_relaxer.process(prot=outs[key]["unrelaxed_protein"])
      pred_output_path = os.path.join(output_dir,f'{args.prefix}_{prefix}_relaxed.pdb')
      with open(pred_output_path, 'w') as f:
        f.write(relaxed_pdb_lines)

############################################################
print(f"model rank based on {rank_by}")
for n,key in enumerate(model_rank):
  print(f"rank_{n+1}_{key} {rank_by}:{outs[key][rank_by]:.2f}")
  if use_ptm and save_pae_json:
    pae = outs[key]["pae"]
    max_pae = pae.max()
    # Save pLDDT and predicted aligned error (if it exists)
    pae_output_path = os.path.join(output_dir,f'rank_{n+1}_{key}_pae.json')
    # Save predicted aligned error in the same format as the AF EMBL DB
    rounded_errors = np.round(np.asarray(pae), decimals=1)
    indices = np.indices((len(rounded_errors), len(rounded_errors))) + 1
    indices_1 = indices[0].flatten().tolist()
    indices_2 = indices[1].flatten().tolist()
    pae_data = json.dumps([{
        'residue1': indices_1,
        'residue2': indices_2,
        'distance': rounded_errors.flatten().tolist(),
        'max_predicted_aligned_error': max_pae.item()
    }],
                          indent=None,
                          separators=(',', ':'))
    with open(pae_output_path, 'w') as f:
      f.write(pae_data)

#@title Extra plots
if use_ptm:
  print("predicted alignment error")
  cf.plot_paes([outs[k]["pae"] for k in model_rank],dpi=dpi)
  plt.savefig(os.path.join(output_dir,f'predicted_alignment_error.png'), bbox_inches = 'tight', dpi=np.maximum(200,dpi))
  plt.show()

print("predicted contacts")
cf.plot_adjs([outs[k]["adj"] for k in model_rank],dpi=dpi)
plt.savefig(os.path.join(output_dir,f'predicted_contacts.png'), bbox_inches = 'tight', dpi=np.maximum(200,dpi))
plt.show()

print("predicted distogram")
cf.plot_dists([outs[k]["dists"] for k in model_rank],dpi=dpi)
plt.savefig(os.path.join(output_dir,f'predicted_distogram.png'), bbox_inches = 'tight', dpi=np.maximum(200,dpi))
plt.show()

print("predicted LDDT")
cf.plot_plddts([outs[k]["plddt"] for k in model_rank], Ls=Ls, dpi=dpi)
plt.savefig(os.path.join(output_dir,f'predicted_LDDT.png'), bbox_inches = 'tight', dpi=np.maximum(200,dpi))
plt.show()

# add settings file
settings_path = os.path.join(output_dir,"settings.txt")
with open(settings_path, "w") as text_file:
  text_file.write(f"sequence={ori_sequence}\n")
  text_file.write(f"msa_method={msa_method}\n")
  text_file.write(f"homooligomer={homooligomer}\n")
  text_file.write(f"pair_msa={pair_msa}\n")
  text_file.write(f"max_msa={max_msa}\n")
  text_file.write(f"cov={cov}\n")
  text_file.write(f"use_amber_relax={use_amber_relax}\n")
  text_file.write(f"use_turbo={use_turbo}\n")
  text_file.write(f"use_ptm={use_ptm}\n")
  text_file.write(f"rank_by={rank_by}\n")
  text_file.write(f"num_models={num_models}\n")
  text_file.write(f"subsample_msa={subsample_msa}\n")
  text_file.write(f"num_samples={num_samples}\n")
  text_file.write(f"num_ensemble={num_ensemble}\n")
  text_file.write(f"max_recycles={max_recycles}\n")
  text_file.write(f"tol={tol}\n")
  text_file.write(f"is_training={is_training}\n")
  text_file.write(f"use_templates=False\n")
  text_file.write(f"-------------------------------------------------\n")

  for n,key in enumerate(model_rank):
    line = f"rank_{n+1}_{key} pLDDT:{outs[key]['pLDDT']:.2f}" + f" pTMscore:{outs[key]['pTMscore']:.4f}" if use_ptm else ""
    text_file.write(line+"\n")
