from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pprint import pprint
import random
from glob import glob
from collections import defaultdict, Counter
import re

def load(infile, format='fasta'):
    ''' Takes path to input file and file format (default = fasta);
    returns list of biopython sequence objects'''

    f = open(infile, 'r')
    try:
        seqs = [i for i in SeqIO.parse(f, format)]
    except ValueError: ## If no format or wrong format provided, try and infer from file extension
        extension = infile.split('.')[-1]
        common_formats = {'phy': 'phylip', 'phyx': 'phylip', 'gb': 'genbank', 'fst': 'fasta', 'fasta': 'fasta'}
        assert extension in common_formats, 'ERROR: provide file format'
        seqs = [i for i in SeqIO.parse(f, common_formats[extension])]

    f.close()

#     print 'Loaded %d sequences from file %s'%(len(seqs), infile)
    return seqs

def load_reference(reference_file):
    reference_seq = SeqIO.read(reference_file, 'genbank')
    gene_locs = {}

    genome_annotation = reference_seq.features
    for f in reference_seq.features:
        if 'gene' in f.qualifiers:
            gene_locs[f.qualifiers['gene'][0]] = (int(f.location.start), int(f.location.end))#FeatureLocation(start=f.location.start, end=f.location.end, strand=1)
    return gene_locs, reference_seq

def categorize(seqlist, fn=None, fields=[1], sep='|'):
    '''[sequenceobjects] --> { 'category1': [sequenceobjects in category1]}'''

    if fn == None:
        fn = lambda x: tuple([x.split(sep)[i] for i in fields]) ## Inefficient, should be changed later

    categories = defaultdict(list)
    for i in seqlist:
        categories[fn(i.description)].append(i)
    return dict(categories)

def find_duplicates(seqs):
    '''Takes a list of aligned sequence objects. Checks for pairs of sequences with identical sequences at all non-gapped sites'''

    assert len(set( [len(i.seq) for i in seqs] )) == 1, "ERROR: Looks for duplicates in aligned sequences. Must be the same length."

    duplicates = [] ## Is there a better way to match groups of sequences that are identical?

    for i in range(len(seqs)): # compare each pair of sequences only once
        for j in range(i+1, len(seqs)):
            zipseqs = zip(str(seqs[i].seq), str(seqs[j].seq))
            zipseqs = [(iseq, jseq) for (iseq, jseq) in zipseqs if iseq != '-' and jseq != '-'] #restrict to nongapped sites
            identicalsites = [(iseq, jseq) for (iseq, jseq) in zipseqs if iseq == jseq]
            if len(identicalsites) == len(zipseqs):
                duplicates.append((seqs[i],seqs[j]))

    if duplicates == []:
        return False
    else:
        return duplicates


    ## Alt. ideas:
    ## Just return a boolean using np.unique; this is nonideal because it won't tell you WHICH sequences are duplicates
    ## Use set logic to test number of unique sequences present; this is nonideal because it doesn't account for identical-except-for-gaps sequences.

def summarize_seqs(seqs):
    ''' Get min, avg, and max values for basic stats about a set of sequences. Currently just nt data, expand later. '''
    stats = {
    'n_sequences': len(seqs),
    'aligned': False,
    'duplicates': None,
    'tiny_seqs': [],
    'unambiguous_length': [],
    'non_actg': [],
    'gc_content': [] }

    if len(set( [len(i.seq) for i in seqs] )) == 1:
        aligned = True
        stats['duplicates'] = find_duplicates(seqs)
    for i in seqs:
        unambiguous_seq = [ nt for nt in str(i.seq).lower() if nt in ['a', 'c', 't', 'g', 'u'] ]
        if len(unambiguous_seq) < 50:
            tiny_seqs.append(i)
        stats['unambiguous_length'].append(len(unambiguous_seq))
        stats['non_actg'].append((float(len(i.seq) - len(unambiguous_seq)) / float(len(i.seq))))
        stats['gc_content'].append(sum([1. if (nt=='g' or nt=='c') else 0. for nt in unambiguous_seq ]) / float(len(unambiguous_seq)))

    def summarize_vals(vals):
        avg = float(sum(vals)) / float(len(vals))
        return ( min(vals), avg, max(vals) )

    for key in ['unambiguous_length', 'non_actg', 'gc_content']:
        stats[key] = summarize_vals(stats[key])

    pprint(stats)
    return stats

def validate(seqs, halt=False):
    '''Screens for common problems in nucleotide sequence sets.
    Raises assertion error if any of the following detected:
    * short sequences < 50nt long
    * >30% ambiguous sites
    * wide range of GC content
    For alignments, also screens for duplicates.
    Will eventually also screen for recombination, abnormal divergence, and back translation.
    '''
    stats = summarize_seqs(seqs)
    problems = { k: False for k in stats.keys() }

    while True not in problems.values():
        problems['duplicates'] = stats['duplicates'] # no duplicates
        problems['tiny_seqs'] = (stats['tiny_seqs'] != []) # no tiny sequences
        problems['non_actg'] = (stats['non_actg'][-1] >= 0.3) # max 30% ambiguous sites

        # max range of 10% difference between most/least-GC rich sequences in set
        problems['gc_content'] = (stats['gc_content'][-1] - stats['gc_content'][0] >= 0.1)
#         if stats['aligned'] == True: # If no sequence set problems, look for alignment problems if relevant
            # @ sites where both valid bases == --> gytis's r-squared code @ flu b repo
            # look for too many mutations from consensus
            # back translation: is each aa encoded by >= 1 aa

        print 'No problems found in sequence set. Screened for duplicates, tiny sequences, ambiguous sites, and GC content. If aligned, also screened for duplicates.'
        break
    else:
        assert True not in problems.values(), ('Found this red flag in your sequence set:', problems)

def count(seqs, fn=None, fields=[1], sep='|'):
    '''
    Return dictionary of counts of sequences for each category. Categories pulled from description via either:
    1 - Passed lambda function. if provided, overrides fields and sep arguments.
    OR
    2 - (Combination of) value(s) of field(s) split by sep
    '''

    if fn == None:
        fn = lambda x: tuple([x.split(sep)[i] for i in fields]) ## Inefficient, should be changed later

    counts = Counter([fn(i.description) for i in seqs])
    return dict(counts) # return counts as a normal dictionary

def rename(seqs, fn=None, fields=[0], sep='|'):
    '''
    Takes and returns a list of biopython sequence objects.
    Transforms sequence descriptions and names IN PLACE according to either:
    1 - Passed lambda function. If provided, overrides fields and sep arguments.
    OR
    2 - List of sep-delimited field indices to keep.
    E.g., if fields=[0,2], sep='|', then >accession|country|date yields >accession|date
    '''
    if fn == None:
        fn = lambda x: sep.join([x.split(sep)[i] for i in fields]) ## Inefficient, should be changed later
    for i in seqs:
        i.description = fn(i.description)
        i.name = fn(i.name)
        i.id = fn(i.id) ## Double check this still works the same way for non-fasta files
    return seqs

def subsample(seqs, k, fn = None, fields=[1], sep='|', multiset_method = 'hierarchy'):
    '''
    seqs: either [biopython sequence objects], or [ [biopython sequence objects, set1], [more sequence objects, set2], ...]
    k: number of sequences to sample from each category. Categories defined based on sequence description via either:
    1 - Passed lambda function. if provided, overrides fields and sep arguments. ## add date parsing later
    OR
    2 - (Combination of) value(s) of field(s) split by sep

    If seqs is a list of lists, handle multiple datasets as specified in multiset_method by either:
    hierarchy: preferentially pull from the first sequence list, then the second, ... until k is reached for each category.
    pool: pool all the sequence lists together, sample randomly
    '''
    if fn == None:
        fn = lambda x: tuple([x.split(sep)[i] for i in fields]) ## Inefficient, should be changed later

    def sample(seqlist, quotas):
        '''
        seqlist: [sequence objects], quotas: {'category1': 20, ...} -->
        {'category1': [20 sequences from category1 or all available sequences from category1 if < 20 available]}
        '''
        categorized = categorize(seqlist, fn=fn)
        the_chosen_ones = {}

        for category, seqsubset in categorized.items():
            quota = quotas[category]

            if quota > 0:
                the_chosen_ones[category] = random.sample(seqsubset, min(quota, len(seqsubset)))
        return the_chosen_ones

    if multiset_method == 'pool' and not isinstance(seqs[0], SeqRecord): # Pool sequence lists together if indicated
        pooled_sequences = []
        for seqlist in seqs:
            pooled_sequences.extend(seqlist)
        seqs = pooled_sequences

    if isinstance(seqs[0], SeqRecord): # If only handed one set of sequences to sample from (or they've been pooled)
        all_categories = set([fn(i.description) for i in seqs])
        quotas = { category: k for category in all_categories } # just take however many we have if <= k
        return sample(seqs, quotas)

    else: # Hierarchical sampling
        all_categories = []
        for seqlist in seqs:
            all_categories.extend([fn(i.description) for i in seqlist])
        quotas = { category: k for category in set(all_categories) }

        subsampled_seqs = Counter()
        for seqlist in seqs:
            subsampled_seqs += Counter(sample(seqlist, dict(quotas))) # i.e., extend for {hash: list} dictionaries
            quotas = Counter(quotas) - subsampled_seqs
        return dict(subsampled_seqs)

				def convert_coordinate(refseq, compareseq, coordinate):
				    '''
				    In: ungapped reference sequence, gapped (aligned) compare sequence, position to convert
				    Out: coordinate of corresponding condition in the aligned sequence
				    '''
				    coordinate = coordinate - 1 # Adjust for python coordinates
				    if isinstance(refseq, SeqRecord):
				        refseq = str(refseq.seq)
				    if isinstance(compareseq, SeqRecord):
				        compareseq = str(compareseq.seq)
				    #check to make sure we have at least 100bp of downstream sequence, no more than 40bp of which are gaps, to match on
				    assert len(refseq[coordinate:]) >= 100 and refseq[coordinate:coordinate+100].count('-')<40, 'ERROR: Not enough downstream context available for '+str(coordinate+1)
				    #check to make sure the given coordinate doesn't correspond to a gap in the reference sequence
				    assert refseq[coordinate] != '-', 'ERROR! Coordinate '+str(coordinate+1)+' is a gap in the reference sequence.'

				    reference = refseq[coordinate:coordinate+100].replace('-', '')
				    refpattern = '-*'.join(list(reference)) # Match in new sequence while ignoring gaps
				    matchlist = re.findall(refpattern, compareseq, flags=re.IGNORECASE) # Check to make sure we get one and only one match
				    assert len(matchlist) == 1, 'ERROR: found %d matches for coordinate %d'%(len(matchlist), coordinate+1)
				    return re.search(refpattern, compareseq, flags=re.IGNORECASE).start()+1 #return the converted coordinate (adjusted back to genomic coordinates)

				def split_alignment(alignment, positions, flanking = True, prune=0.25, write = True, ofilestem='split'):
				    '''
				    In: alignment path or list of seq objects, integer position to split the alignment along, or a list of them.

				    Args:
				    flanking: if True, includes the segments [0, firstposition], [lastposition, end].
				    prune: if float between 0 and 1, removes any sequences with <= that proportion of non-gapped/N sites
				    ofilestem, write: if True, write each segmet to cwd$ ofilestem_start_end.fasta

				    N.B.: Splitting is python-esque (includes start point, does not include endpoint)

				    Out: Returns dictionary of {(start, end): [split alignment objects]};
				    optionally writes each segment to ofilestem_start_end.fasta
				    '''
				    ## Handle input data types
				    if isinstance(alignment, str):
				        alignment = load(alignment)
				    else:
				        assert isinstance(alignment, list) and all([isinstance(i, SeqRecord) for i in alignment]) # If not a path name, should be a list of Seq objects
				    assert len(set([len(i.seq) for i in alignment])) == 1, 'ERROR: All sequences must be the same length (aligned)' ## All sequences are the same length

				    if isinstance(positions, int):
				        positions = [positions]
				        assert flanking == True, 'ERROR: only one position given. If you want [0, position], [position, end], then run with flanking=True'
				    else:
				        assert isinstance(positions, list) and all([isinstance(i, int) for i in positions])

				    positions.sort()

				    ## Add endpoints if flanking=True
				    if flanking:
				        positions.insert(0, 0)
				        positions.append(len(alignment[0]))

				    ## Split alignment
				    splits = {}
				    for i, p in enumerate(positions[:-1]):
				        start = positions[i]
				        end = positions[i+1] ## like python splitting, does not include endpoint

				        if all([isinstance(prune, float), prune >= 0., prune <= 1.]): # if removing uninformative sequences....
				            splits['%d_%d'%(start, end)] = []
				            for s in alignment:
				                ambiguity_codes = []
				                prop_ambiguous = sum([ 1. if site in ambiguity_codes else 0. for site in str(s.seq)]) / float(len(s.seq))
				                if prop_ambiguous <= prune:
				                    splits['%d_%d'%(start, end)].append(s)
				        else:
				            splits['%d_%d'%(start, end)] = [s[start:end] for s in alignment]

				    ## Optionally write out sequences
				    if write == True:
				        for ((start, end), seqs) in splits.iteritems():
				            SeqIO.write(seqs, '%s_%s_%s.fasta'%(ofilestem, start, end), 'fasta')

				    return splits


				def extract_gene(alignment, genes, reference_file=None, gene_locs=None, reference_seq=None, prune=0.25, write = True):
				    '''
				    In: alignment path or list of seq objects, name(s) of genes to extract (str or list),
				    reference file path (overrides) or output from load_reference,
				    max ambiguous proportion for pruning in split_alignment, boolean for whether to write each gene alignment to file.

				    Out: Returns dictionary of {gene: [split alignment objects]};
				    optionally writes each segment to gene_start_end.fasta
				    '''

				    ## Handle input data types
				    if isinstance(alignment, str):
				        alignment = load(alignment)
				    else:
				        assert isinstance(alignment, list) and all([isinstance(i, SeqRecord) for i in alignment]) # If not a path name, should be a list of Seq objects
				    assert len(set([len(i.seq) for i in alignment])) == 1, 'ERROR: All sequences must be the same length (aligned)' ## All sequences are the same length

				    if isinstance(genes, str):
				        genes = [ genes ]

				    if reference_file == None:
				        assert isinstance(gene_locs, dict) and isinstance(reference_seq, SeqRecord), "ERROR: Either provide path to reference file or output from load_reference"
				    else:
				        gene_locs, reference_seq = load_reference(reference_file)
				    assert reference_seq.name in [i.name for i in alignment], "ERROR: Reference seq must be in alignment"
				    aligned_ref_seq = [i for i in alignment if i.name == reference_seq.name][0]

				    extracted_genes = {}
				    for gene in genes:
				        assert gene in gene_locs, ("ERROR: gene names must be in reference file. Try: ", gene_locs.keys())

				        start, end = [convert_coordinate(reference_seq, aligned_ref_seq, i) for i in gene_locs[gene]]
				        extracted_genes[gene] = split_alignment(alignment, [start, end], flanking=False, prune=prune, write=write, ofilestem=gene)['%s_%s'%(start, end)]

				    return extracted_genes
