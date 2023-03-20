import numpy as np
import pandas as pd
import pyfaidx
import pyranges


def gene_landscape(df):
    """
    inspired by: https://www.nature.com/articles/nature16994 Extended Data Figure 5a
    
    note: this implicitely assumes equal stepsize between timepoints
    """
    new_df = pd.DataFrame(index=df.index, columns=["late", "mid-development"])
    
    timepoints = df.columns
    
    late = pd.Series(np.linspace(0, 1, len(timepoints)), index=timepoints)
    if len(timepoints) % 2 == 0:
        switch = pd.Series(
            np.concatenate([
                np.linspace(0, 1, len(timepoints) // 2),
                np.linspace(1, 0, len(timepoints) // 2)]), 
            index=timepoints)
    else:
        switch = pd.Series(
            np.concatenate([
                np.linspace(0, 1, (len(timepoints) // 2) + 1),
                np.linspace(1, 0, len(timepoints) // 2)]), 
            index=timepoints)
    
    new_df["late"] = df.corrwith(late, axis=1)
    new_df["mid-development"] = df.corrwith(switch, axis=1)
    
    return new_df.dropna()


def pyranges2psl(pyranges_peaks, psl, genome):
    """
    Converts a dictionary of peaks to a psl file. The keys of the dictionary are ignored,
    but the values are pyranges which are then converted to the psl format.
    """
    # we read the whole genome into memory as pyfaidx doesn't work well with gzipped genomes
    genome = pyfaidx.Fasta(genome)
    genome = {contig: seq.__str__() for contig, seq in genome.items()}

    with open(psl, "w") as f:
        for timepoint, enhancerss in pyranges_peaks.items():
            for _, enhancers in enhancerss:
                for _, row in enhancers.iterrows():
                    sequence = genome[row.Chromosome][int(row.Start):int(row.End)]
                    contig = row.Chromosome
                    start = row.Start
                    stop = row.End

                    res = [
                        len(sequence), 
                        0, 
                        0, 
                        sequence.lower().count("n"), 
                        0, 
                        0, 
                        0, 
                        0, 
                        "+", 
                        f"{contig}:{start}-{stop}", 
                        len(sequence), 
                        0, 
                        len(sequence),
                        contig,
                        len(genome[contig]),
                        start,
                        stop,
                        1,
                        f"{len(sequence)},",
                        f"0,",
                        f"{start},",
                    ]
                    res = [str(x) for x in res]
                    f.write("\t".join(res) + "\n")


def jaccard_peaks(dmel, dvir, peak_conversion, nr_peaks=None, timepoints=["TP1", "TP2", "TP3", "TP4", "TP5"]):
    """
    Function to get the jaccard index of overlap of two peaksets (dmel & dvir)
    based on orthology (peak_conversion).
    """
    res = {}
    for tp in timepoints:
        # load the peaks
        tp_specific_dmel_sample = dmel[tp]
        tp_specific_dvir_sample = dvir[tp]
        if nr_peaks is not None:
            tp_specific_dmel_sample = dmel[tp].sample(nr_peaks)
            tp_specific_dvir_sample = dvir[tp].sample(nr_peaks)

        # look up the virilis peaks in the conversion table
        tp_specific_dvir_sample_set = set()
        for chrom, locs in tp_specific_dvir_sample:
            if isinstance(chrom, tuple):
                chrom = chrom[0]
            for _, row in locs.iterrows():
                tp_specific_dvir_sample_set.add(f"{chrom}:{row.Start}-{row.End}")
        _peak_conversion = peak_conversion[peak_conversion[9].isin(tp_specific_dvir_sample_set)]

        # prepare the converted virilis peaks for pyranges
        peaks_dvir = {"Chromosome": [], "Start": [], "End": [], "origin": []}
        for i, row in _peak_conversion.iterrows():
            peaks_dvir["Chromosome"].append(row.loc[13])
            peaks_dvir["Start"].append(row.loc[15])
            peaks_dvir["End"].append(row.loc[16])
            peaks_dvir["origin"].append(row.loc[9])

        # load the converted peaks into pyranges
        peaks_dvir = pyranges.from_dict(peaks_dvir)

        # get the intersection between melanogaster peaks and converted virilis peaks
        conserved_enhancer_peaks = tp_specific_dmel_sample.intersect(peaks_dvir)
        n_dmel = len(tp_specific_dmel_sample)
        n_dvir = len(tp_specific_dvir_sample)
        intersect = len(conserved_enhancer_peaks)
        res[tp] = (intersect) / (n_dmel + n_dvir - intersect)
    
    return res



def get_peaks(bedfiles, timepoint_specific=False):
    """
    
    """
    dmel_peaks = {}
    dvir_peaks = {}
    for i, tp in enumerate(["TP1", "TP2", "TP3", "TP4", "TP5"]):
        peaks_dmel = {"Chromosome": [], "Start": [], "End": []}
        peaks_dvir = {"Chromosome": [], "Start": [], "End": []}
        # 
        if timepoint_specific:
            for dmel, dvir in zip(bedfiles[0][:i] + bedfiles[0][i+1:],
                                  bedfiles[1][:i] + bedfiles[1][i+1:]):
                dmel = pd.read_table(dmel, usecols=[0, 1, 2], names=["chrom", "start", "stop"])
                dvir = pd.read_table(dvir, usecols=[0, 1, 2], names=["chrom", "start", "stop"])
                for _, row in dmel.iterrows():
                    peaks_dmel["Chromosome"].append(row.chrom)
                    peaks_dmel["Start"].append(row.start)
                    peaks_dmel["End"].append(row.stop)
                for _, row in dvir.iterrows():
                    peaks_dvir["Chromosome"].append(row.chrom)
                    peaks_dvir["Start"].append(row.start)
                    peaks_dvir["End"].append(row.stop)

            # then we read the dictionary into pyranges and merge all peaks
            tp_specific_dmel = pyranges.from_dict(peaks_dmel).merge()
            tp_specific_dvir = pyranges.from_dict(peaks_dvir).merge()

            # then we read the time point peaks and find the nearest
            tp_specific_dmel = pyranges.PyRanges(df=pd.read_table(bedfiles[0][i], usecols=range(3), index_col=False, names=["Chromosome", "Start", "End"])).nearest(tp_specific_dmel)
            tp_specific_dvir = pyranges.PyRanges(df=pd.read_table(bedfiles[1][i], usecols=range(3), index_col=False, names=["Chromosome", "Start", "End"])).nearest(tp_specific_dvir)

            # and only keep the time point peaks which are not overlapping (distance > 0)
            dmel_peaks[tp] = pyranges.PyRanges(df=tp_specific_dmel.df[tp_specific_dmel.df["Distance"] > 0])
            dvir_peaks[tp] = pyranges.PyRanges(df=tp_specific_dvir.df[tp_specific_dvir.df["Distance"] > 0])
        #
        else:
            dmel = pd.read_table(bedfiles[0][i], usecols=[0, 1, 2], names=["chrom", "start", "stop"])
            dvir = pd.read_table(bedfiles[1][i], usecols=[0, 1, 2], names=["chrom", "start", "stop"])
            for _, row in dmel.iterrows():
                peaks_dmel["Chromosome"].append(row.chrom)
                peaks_dmel["Start"].append(row.start)
                peaks_dmel["End"].append(row.stop)
            for _, row in dvir.iterrows():
                peaks_dvir["Chromosome"].append(row.chrom)
                peaks_dvir["Start"].append(row.start)
                peaks_dvir["End"].append(row.stop)

            # then we read the dictionary into pyranges and merge all peaks
            dmel_peaks[tp] = pyranges.from_dict(peaks_dmel).merge()
            dvir_peaks[tp] = pyranges.from_dict(peaks_dvir).merge()
    return dmel_peaks, dvir_peaks


def jensen_shannon_distance(x, y):
    # kullback leibler divergence
    kl_div = lambda x,y: np.nansum(x * np.log2(x/y))

    # normalize to sum to one
    x, y = x / x.sum(), y / y.sum()

    m = 0.5 * (x + y)
    return np.sqrt(0.5 * (kl_div(x, m) + kl_div(y, m)))


def gen_orthodb(orthofile, db=":memory:"):
    """
    generate a ortholog database
    
    TODO
    """
    import os
    import sqlite3
    import urllib.request, json 

    orthogroups = pd.read_csv(orthofile, sep="\t", comment="#", index_col=0)
    if os.path.exists(db):
        os.remove(db)
    conn = sqlite3.connect(db)

    # make the orthogroup table
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS orthogroups
        (
        id integer PRIMARY KEY, 
        orthogroup text UNIQUE NOT NULL
        )
        """
    )

    # make the assemblies table
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS assemblies
        (
        id integer PRIMARY KEY, 
        assembly text
        )
        """
    )

    # make the genes table
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS genes
        (
        id integer PRIMARY KEY,
        gene_name text,
        gene_id text,
        assembly NOT NULL,
        orthogroup,
        FOREIGN KEY (orthogroup) REFERENCES orthogroups (id),
        FOREIGN KEY (assembly) REFERENCES assemblies (assembly)
        )
        """
    )
    
    # add all orthogroups to the table
    for orthogroup in orthogroups.index:
        conn.execute(f"INSERT INTO orthogroups VALUES(NULL, '{orthogroup}')")
    conn.execute(f"INSERT INTO orthogroups VALUES(NULL, 'UNASSIGNED')")

    # add all assemblies
    all_assemblies = [assembly.rstrip(".pep") for assembly in orthogroups.columns if assembly.endswith(".pep")]
    for assembly in all_assemblies:
        conn.execute(f"INSERT INTO assemblies VALUES(NULL, '{assembly}')")

    # finally add all genes per species with a connection to the orthogroup
    for assembly in all_assemblies:
        for orthogroup, genes in zip(
            orthogroups.index, orthogroups[f"{assembly}.pep"]
        ):
            if isinstance(genes, float) and np.isnan(genes):
                continue

            for gene in genes.split(", "):
                gene_name, gene_id = gene.split("|")
                if gene_name != ".":
                    conn.execute(
                        f"INSERT INTO genes VALUES(NULL, '{gene_name}', '{gene_id}', '{assembly}', '{orthogroup}')"
                    )
                else:
                    conn.execute(
                        f"INSERT INTO genes VALUES(NULL, NULL, '{gene_id}', '{assembly}', '{orthogroup}')"
                    )

    conn.commit()
    cur = conn.cursor()

    cur.execute("CREATE INDEX idx_name ON genes (gene_name, assembly)")
    cur.execute("CREATE INDEX idx_id ON genes (gene_id, assembly)")
    cur.execute("CREATE INDEX idx_gene_name ON genes (gene_name)")
    cur.execute("CREATE INDEX idx_gene_id ON genes (gene_id)")
    cur.execute("CREATE INDEX idx_assembly_gene_name ON genes (assembly, gene_name)")
    cur.execute("CREATE INDEX idx_assembly_gene_id ON genes (assembly, gene_id)")
    cur.execute("CREATE INDEX idx_orthogroup ON orthogroups (orthogroup)")
    cur.execute("CREATE INDEX idx_orthogroup_assembly ON genes (orthogroup, assembly)")
    cur.execute("CREATE INDEX idx_assembly_orthogroup ON genes (assembly, orthogroup)")
    conn.commit()
    
    return cur


    
def pairwise_spearman(df1, df2):
    corr_matrix = np.zeros((len(df1.columns), len(df2.columns)))
    
    df1 = df1.rank(axis="rows")
    df2 = df2.rank(axis="rows")
    for k, x_col in enumerate(df1.columns):
        for l, y_col in enumerate(df2.columns):
            corr_matrix[k, l] = np.corrcoef(df1[x_col], df2[y_col])[0, 1]
#     print(corr_matrix)
    return corr_matrix
    
def pairwise_jensen_shannon(df1, df2):
    corr_matrix = np.zeros((len(df1.columns), len(df2.columns)))
    
    for k, x_col in enumerate(df1.columns):
        for l, y_col in enumerate(df2.columns):
            corr_matrix[k, l] = jensen_shannon_distance(df1[x_col], df2[y_col])
    return corr_matrix

def pairwise_pearson(df1, df2):
    corr_matrix = np.zeros((len(df1.columns), len(df2.columns)))
    
    for k, x_col in enumerate(df1.columns):
        for l, y_col in enumerate(df2.columns):
            corr_matrix[k, l] = np.corrcoef(df1[x_col], df2[y_col])[0, 1]
    return corr_matrix
