import os
import re

# Get rid of intel-dependent modules using "module purge" command and
# then capture list of available modules in single column format
stream = os.popen("module purge; module avail -t 2>&1 >/dev/null")

# Print header
print 'Name,Version,Default,Type,Domain,Handle Type,Handle Key,Description'

# Set domain dictionary
domain_dict = {
'amber':               'Chemistry',
'apbs':                'Chemistry',
'atlas':               'Mathematics',
'beagle':              'Phylogenetics',
'beast':               'Phylogenetics',
'beast2':              'Phylogenetics',
'bedtools':            'Bioinformatics',
'bismark':             'Bioinformatics',
'bowtie':              'Bioinformatics',
'bamtools':            'Bioinformatics',
'biopython':           'Bioinformatics',
'blast':               'Bioinformatics',
'blat':                'Bioinformatics',
'bowtie2':             'Bioinformatics',
'bx-python':           'Bioinformatics',
'cufflinks':           'Bioinformatics',
'dendropy':            'Bioinformatics',
'edena':               'Bioinformatics',
'fastqc':              'Bioinformatics',
'fastx':               'Bioinformatics',
'genomeanalysistk':    'Bioinformatics',
'gmap_gsnap':          'Bioinformatics',
'htseq':               'Bioinformatics',
'idba-ud':             'Bioinformatics',
'matt':                'Bioinformatics',
'mirdeep2':            'Bioinformatics',
'miso':                'Bioinformatics',
'picard':              'Bioinformatics',
'plink':               'Bioinformatics',
'pysam':               'Bioinformatics',
'randfold':            'Bioinformatics',
'samtools':            'Bioinformatics',
'soapdenovo':          'Bioinformatics',
'soapsnp':             'Bioinformatics',
'spades':              'Bioinformatics',
'squid':               'Bioinformatics',
'tophat':              'Bioinformatics',
'trimmomatic':         'Bioinformatics',
'trinity':             'Bioinformatics',
'velvet':              'Bioinformatics',
'viennarna':           'Bioinformatics',
'boost':               'General',
'cp2k':                'Chemistry/Bio',
'cpmd':                'Molecular dynamics',
'ddt':                 'Debugging/Profiling',
'dppdiv':              'Bioinformatics',
'fasttree':            'Bioinformatics',
'fsa':                 'Bioinformatics',
'gamess':              'Quantum Chemistry',
'garli':               'Bioinformatics',
'gaussian':            'Chemistry',
'gromacs':             'Molecular dynamics',
'hadoop':              'Data analysis',
'idl':                 'Data analysis/Visualization',
'lammps':              'Molecular dynamics',
'mafft':               'Bioinformatics',
'matlab':              'Mathematics',
'mopac':               'Chemistry',
'mrbayes':             'Bioinformatics',
'namd':                'Molecular dynamics',
'nwchem':              'Chemistry',
'octave':              'General',
'phylobayes':          'Bioinformatics',
'q-chem':              'Chemistry',
'quantum espresso':    'Chemistry',
'r':                   'Statistical Computing',
'raxml':               'Bioinformatics',
'scipy':               'General',
'spark':               'Data analysis',
'visit':               'Visualization',
'weka':                'Data Analytics'
}

# Set type dictionary
type_dict = {
'amber':               'package', 
'apbs':                'package', 
'atlas':               'library', 
'beagle':              'library', 
'beast':               'package', 
'beast2':              'package', 
'bedtools':            'package', 
'bismark':             'package', 
'bowtie':              'package', 
'bamtools':            'toolkit', 
'biopython':           'package', 
'blast':               'package', 
'blat':                'package', 
'bowtie2':             'package', 
'bx-python':           'library', 
'cufflinks':           'package', 
'dendropy':            'library', 
'edena':               'package', 
'fastqc':              'tool',    
'fastx':               'toolkit', 
'genomeanalysistk':    'toolkit', 
'gmap_gsnap':          'package', 
'htseq':               'package', 
'idba-ud':             'package', 
'matt':                'package', 
'mirdeep2':            'tool',    
'miso':                'package', 
'picard':              'toolkit', 
'plink':               'toolkit', 
'pysam':               'module',  
'randfold':            'package', 
'samtools':            'toolkit', 
'soapdenovo':          'package', 
'soapsnp':             'tool',    
'spades':              'tool',    
'squid':               'package', 
'tophat':              'program', 
'trimmomatic':         'tool',    
'trinity':             'program', 
'velvet':              'package', 
'viennarna':           'package', 
'boost':               'library', 
'cp2k':                'program', 
'cpmd':                'program', 
'ddt':                 'package', 
'dppdiv':              'package', 
'fasttree':            'package', 
'fsa':                 'package', 
'gamess':              'package', 
'garli':               'program', 
'gaussian':            'package', 
'gromacs':             'package', 
'hadoop':              'library', 
'idl':                 'package', 
'lammps':              'package', 
'mafft':               'program', 
'matlab':              'program', 
'mopac':               'program', 
'mrbayes':             'program', 
'namd':                'program', 
'nwchem':              'package', 
'octave':              'package', 
'phylobayes':          'package', 
'q-chem':              'package', 
'quantum espresso':    'package', 
'r':                   'platform',
'raxml':               'tool',    
'scipy':               'library', 
'spark':               'library', 
'visit':               'package', 
'weka':                'package' 
}

# Set handle_type dictionary
handle_type_dict = {
'amber':               'module', 
'apbs':                'module', 
'atlas':               'module', 
'beagle':              'module', 
'beast':               'module', 
'beast2':              'module', 
'bedtools':            'module', 
'bismark':             'module', 
'bowtie':              'module', 
'bamtools':            'module', 
'biopython':           'module', 
'blast':               'module', 
'blat':                'module', 
'bowtie2':             'module', 
'bx-python':           'module', 
'cufflinks':           'module', 
'dendropy':            'module', 
'edena':               'module', 
'fastqc':              'module',    
'fastx':               'module', 
'genomeanalysistk':    'module', 
'gmap_gsnap':          'module', 
'htseq':               'module', 
'idba-ud':             'module', 
'matt':                'module', 
'mirdeep2':            'module',    
'miso':                'module', 
'picard':              'module', 
'plink':               'module', 
'pysam':               'module',  
'randfold':            'module', 
'samtools':            'module', 
'soapdenovo':          'module', 
'soapsnp':             'module',    
'spades':              'module',    
'squid':               'module', 
'tophat':              'module', 
'trimmomatic':         'module',    
'trinity':             'module', 
'velvet':              'module', 
'viennarna':           'module', 
'boost':               'module', 
'cp2k':                'module', 
'cpmd':                'module', 
'ddt':                 'module', 
'dppdiv':              'module', 
'fasttree':            'module', 
'fsa':                 'module', 
'gamess':              'module', 
'garli':               'module', 
'gaussian':            'module', 
'gromacs':             'module', 
'hadoop':              'module', 
'idl':                 'module', 
'lammps':              'module', 
'mafft':               'module', 
'matlab':              'module', 
'mopac':               'module', 
'mrbayes':             'module', 
'namd':                'module', 
'nwchem':              'module', 
'octave':              'module', 
'phylobayes':          'module', 
'q-chem':              'module', 
'quantum espresso':    'module', 
'r':                   'module',
'raxml':               'module',    
'scipy':               'module', 
'spark':               'module', 
'visit':               'module', 
'weka':                'module' 
}

# Set description dictionary
description_dict = {
'amber':               'Package of molecular simulation programs',

'apbs':                'Adaptive Poisson-Boltzmann Software package. Used to analyze the solvation properties of small and macro-molecules.',

'atlas':               'Linear Algebra library.',

'beagle':              'Library for evaluating phylogenetic likelihoods',

'beast':               'Program for Bayesian analysis of molecular sequences',

'beast2':              'Program for Bayesian analysis of molecular sequences',

'bedtools':            'Tools for genomic analysis tasks',

'bismark':             'A tool to map bisulfite converted sequence reads and determine cytosine methylation states',

'bowtie':              'An ultrafast,  memory efficient short read aligner',

'bamtools':            'Fast flexible C++ API and toolkit for reading, writing, and manipulating BAM files',

'biopython':           'Python based tools for biological computation',

'blast':               'Basic local alignment search tool',

'blat':                'Pairwise sequence alignment algorithm',

'bowtie2':             'An ultrafast,  memory efficient short read aligner',

'bx-python':           'python library for rapid implementation of genome scale analyses',

'cufflinks':           'Transcript assembly, differential expression, and differential regulation for RNA-Seq',

'dendropy':            'Python library for phylogenetic computing',

'edena':               'De novo short reads assembler',

'fastqc':              'A quality control tool for high throughput sequence data',

'fastx':               'Collection of command line tools for Short Reads FASTA/FASTQ files preprocessing',

'genomeanalysistk':    'Genome Analysis Toolkit',

'gmap_gsnap':          'Genomic Mapping and Alignment Program for mRNA and EST Sequences',

'htseq':               'Python package that provides infrastructure to process data from high-throughput sequencing assays.',

'idba-ud':             'Iterative De Bruijn Graph De Novo Assembler for Short Reads Sequencing data with Highly Uneven Sequencing Depth',

'matt':                'TBD',

'mirdeep2':            'Tool to discover known and novel miRNAs from deep sequencing data',

'miso':                'Mixture of Isoforms model (MISO) for isoform quantitation using RNA-Seq',

'picard':              'Command line tools for working with next generation sequencing data in the BAM format',

'plink':               'Whole genome association analysis toolset.',

'pysam':               'Python module to read SAM files',

'randfold':            'Minimum free energy of folding randomization test software',

'samtools':            'SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.',

'soapdenovo':          'SOAPdenovo is a novel short-read assembly method that can build a de novo draft assembly for the human-sized genomes. The program is specially de signed to assemble Illumina GA short reads.',

'soapsnp':             'A resequencing utility that can assemble consensus sequence for the genome of a newly sequenced individual based on the alignment of the raw sequencing reads on the known reference.',

'spades':              'A genome assembly algorithm designed for single cell and multi-cells bacterial data sets.',

'squid':               'Bioinformatics grid software',

'tophat':              'TopHat is a fast splice junction mapper for RNA-Seq reads.',

'trimmomatic':         'A flexible read trimming tool for Illumina NGS data',

'trinity':             'A program for the efficient and robust de novo reconstruction of transcriptomes from RNA-seq data.',

'velvet':              'Velvet is an algorithm package that has been designed to deal with de novo genome assembly and short read sequencing alignments.',

'viennarna':           'The ViennaRNA Package consists of a C code library and several stand-alone programs for the prediction and comparison of RNA secondary structures.',

'boost':               'Boost provides free peer-reviewed portable C++ source libraries.',

'cp2k':                'CP2K is a program that performs atomistic and molecular simulations of solid state, liquid, molecular and biological systems.',

'cpmd':                'The CPMD code is a parallelized plane wave / pseudopotentialimplementation of Density Functional Theory, particularly designed for ab-initio molecular dynamics.',

'ddt':                 'Debugging tool to tackle complex multi-threaded or multi-process software problems',

'dppdiv':              'A software tool for estimating species divergence times and lineage-specific substitution rates on a fixed tree topology.',

'fasttree':            'FastTree infers approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences.',

'fsa':                 'FSA is a probabilistic multiple sequence alignment algorithm which uses a distance-based approach to aligning homologous protein, RNA or DNA sequences.',

'gamess':              'The general atomic and molecular electronic structure system (GAMESS) is a general ab initio quantum chemistry package.',

'garli':               'GARLI performs heuristic phylogenetic searches under the General Time Reversible (GTR) model of nucleotide substitution and its submodels, with or without gamma distributed rate heterogeneity and a proportion of invariant sites.',

'gaussian':            'Gaussian provides state-of-the-art capabilities for electronic structure modeling.',

'gromacs':             'GROMACS is a versatile package to perform molecular dynamics',

'hadoop':              'The Apache Hadoop software library is a framework that allows for the distributed processing of large data sets across clusters of computers using simple programming models.',

'idl':                 'IDL is a commercial package that allows users to create high-quality graphic data visualizations, and get statistical, analytical, and numerical information.',

'lammps':              'LAMMPS is a classical molecular dynamics code, and an acronym for Large-scale Atomic/Molecular Massively Parallel Simulator.',

'mafft':               'MAFFT is a multiple sequence alignment program for unix-like operating systems.',

'matlab':              'MATLAB is a high-level language and environment for numerical computation, visualization, and programming. Restricted license, please check with SDSC for details.',

'mopac':               'MOPAC is designed to implement semi-empirical quantum chemistry algorithms.',

'mrbayes':             'MrBayes is a program for Bayesian inference and model choice across a wide range of phylogenetic and evolutionary models.',

'namd':                'NAMD is a parallel molecular dynamics code designed for high-performance simulation of large biomolecular systems.',

'nwchem':              'NWChem provides users with computational chemistry tools that are scalable both in their ability to treat large scientific computational chemistry problems efficiently.',

'octave':              'Octave is a high-level intrepreted languange, primarily intended for numerical computations.',

'phylobayes':          'PhyloBayes is a Bayesian Monte Carlo Markov Chain (MCMC) sampler for phylogenetic reconstruction using protein alignments.',

'q-chem':              'Chemistry software, theoretical chemistry and quantum Chemistry software for research, visualization, quantum calculation and molecular modeling.',

'quantum espresso':    'Integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale.',

'r':                   'R is a free software environment for statistical computing and graphics.',

'raxml':               'RAxML tool for Maximum-likelihood based phylogenetic inference.',

'scipy':               'Open source library of scientific tools.',

'spark':               'Apache Spark is a fast and general engine for large-scale data processing.',

'visit':               'VisIt is an Open Source, interactive, scalable, visualization, animation and analysis tool',

'weka':                'Weka is a collection of machine learning algorithms for data mining tasks.'
}

# Loop over modules
for line in stream:
    default = 'no'
    line = re.sub('\s*', '', line)
    if line == 'dot':
        continue
    if line == 'module-git':
        continue
    if line == 'module-info':
        continue
    if line == 'modules':
        continue
    if line == 'null':
        continue
    if line == 'use.own':
        continue
    if line.find('modulefiles') >= 0:
        continue
    if line.find('(default)') >=0:
        default = 'yes'
        line = re.sub('\(default\)$', '', line)

    handle_key = line
    fields = line.split('/')

    if domain_dict.has_key(fields[0].lower()) :
        mdomain = domain_dict[fields[0].lower()]
    else:
        mdomain = 'undef'
        
    if type_dict.has_key(fields[0].lower()) :
        mtype = type_dict[fields[0].lower()]
    else:
        mtype = 'undef'

    if handle_type_dict.has_key(fields[0].lower()) :
        mhandle_type = handle_type_dict[fields[0].lower()]
    else:
        mhandle_type = 'undef'

    if description_dict.has_key(fields[0].lower()) :
        mdescription = '"' + description_dict[fields[0].lower()] + '"'
    else:
        mdescription = '"' + 'undef' + '"'
        
    print ','.join([fields[0], fields[1], default, mtype, mdomain, mhandle_type, handle_key, mdescription])

