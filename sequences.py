from Bio import SeqIO
from pprint import pprint
import random
from glob import glob
from collections import defaultdict, Counter
import argparse

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
    return seqs

def rename(seqs, fn=None, fields=[0], sep='|'):
    '''
    Takes and returns a list of biopython sequence objects.
    Transforms sequence descriptions and names according to either:
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

def subsample(seqs, k, fn = None, fields=[1], sep='|', multiset_method = 'ordered_preference'):
    '''
    seqs: either a list of biopython sequence objects, or a list of lists of biopython sequence objects
    k: number of sequences to sample from each category

    Sequences sorted into categories based on sequence description via either:
    1 - Passed lambda function. if provided, overrides fields and sep arguments.
    OR
    2 - (Combination of) value(s) of field(s) split by sep

    If seqs is a list of lists, handle multiple datasets as specified in multiset_method by either:
    ordered_preference: preferentially pull from the first sequence list, then the second, ... until k is reached for each category.
    combine: pool all the sequence lists together, sample randomly
    '''
    if fn == None:
        fn = lambda x: tuple([x.split(sep)[i] for i in fields]) ## Inefficient, should be changed later

    def categorize(seqlist):
        categories = defaultdict(list)
        for i in seqlist:
            categories[fn(i.description)].append(i)
        return dict(categories)

    def sample(seqlist, quotas):
        categorized = categorize(seqlist)
        the_chosen_ones = {}

        for category, seqsubset in categorized.items():
            quota = quotas[category]
            if quota > 0:
                the_chosen_ones[category] = random.sample(seqsubset, quota)
        return the_chosen_ones

    if isinstance(seqs[0], Bio.SeqRecord.SeqRecord): # If only handed one set of sequences to sample from
        all_categories = Counter([fn(i.description) for i in seqs])
        quotas = { category: min(k, count) for category, count in all_categories }
        return sample(seqs, quotas)

    else:
        all_categories = []
        for seqlist in seqs:
            all_categories.extend([fn(i.description) for i in seqlist])
        all_categories = Counter(all_categories)
        quotas = { category: min(k, count) for category, count in all_categories }

        for seqlist in seqs:
            
