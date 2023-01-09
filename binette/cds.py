import pyrodigal
import concurrent.futures as cf
import logging
from collections import Counter
import pyfastx
from collections import defaultdict

from tqdm import tqdm

def get_contig_from_cds_name(cds_name):
    return '_'.join(cds_name.split('_')[:-1])

def predict(contigs_iterator, outfaa, threads = 1):
    """Predict open reading frames with Pyrodigal."""

    future_per_contig = {}
    orf_finder = pyrodigal.OrfFinder(meta='meta')

    logging.info(f'Predicting cds sequences with Pyrodigal using {threads} threads.')
    with cf.ProcessPoolExecutor(max_workers=threads) as tpe:
        for seq in contigs_iterator:
            future_per_contig[seq.name] = tpe.submit(orf_finder.find_genes, seq.seq)
    
    contig_to_pyrodigal_genes = {contig_id:future.result() for contig_id, future in  future_per_contig.items()}
    write_faa(outfaa, contig_to_pyrodigal_genes )

    contig_to_genes = {contig_id:[gene.translate() for gene in pyrodigal_genes] for contig_id, pyrodigal_genes in  contig_to_pyrodigal_genes.items()}
    return contig_to_genes


def write_faa(outfaa, contig_to_genes ):
    logging.info('Writting predicted protein sequences.')
    with open(outfaa, 'w') as fl:
        for contig_id, genes in contig_to_genes.items():
            genes.write_translations(fl, contig_id)
        
def parse_faa_file(faa_file):
    contig_to_genes = defaultdict(list)
    for name, seq, _ in pyfastx.Fastx(faa_file):
        contig = get_contig_from_cds_name(name)
        contig_to_genes[contig].append(seq)

    return dict(contig_to_genes)


def get_aa_composition(genes):

    aa_counter = Counter()
    for gene in genes:
        aa_counter += Counter(gene)

    return aa_counter


def get_contig_cds_metadata(contig_to_genes):
    
    contig_to_cds_count = {contig:len(genes) for contig, genes in contig_to_genes.items()}

    contig_to_aa_counter = {contig:get_aa_composition(genes) for contig, genes in tqdm(contig_to_genes.items())}
    logging.info('contig_to_aa_counter done. ')

    contig_to_aa_length = {contig:sum(counter.values()) for contig, counter in contig_to_aa_counter.items()}
    logging.info('contig_to_aa_length done. ')

    return contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length

def get_contig_cds_metadata(contig_to_genes, threads):
    
    contig_to_cds_count = {contig:len(genes) for contig, genes in contig_to_genes.items()}


    contig_to_future = {}
    logging.info(f'Collect contig aa composition using {threads} threads.')
    with cf.ProcessPoolExecutor(max_workers=threads) as tpe:
        for contig, genes in tqdm(contig_to_genes.items()):
            contig_to_future[contig] = tpe.submit(get_aa_composition, genes)
    
    contig_to_aa_counter = {contig:future.result() for contig, future in  tqdm(contig_to_future.items())}
    
    # contig_to_aa_counter = {contig:get_aa_composition(genes) for contig, genes in tqdm(contig_to_genes.items())}
    logging.info('contig_to_aa_counter done. ')
    
    contig_to_aa_length = {contig:sum(counter.values()) for contig, counter in contig_to_aa_counter.items()}
    logging.info('contig_to_aa_length done. ')

    return contig_to_cds_count, contig_to_aa_counter, contig_to_aa_length


